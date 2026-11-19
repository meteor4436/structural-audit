"""Inspect ligand contacts in a PDB file without external dependencies."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

from . import __version__

WATER = {'HOH', 'WAT', 'DOD'}
RADII = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80,
         'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98}


@dataclass(frozen=True)
class Atom:
    record: str
    serial: int
    name: str
    altloc: str
    residue: str
    chain: str
    sequence: str
    insertion: str
    x: float
    y: float
    z: float
    occupancy: float
    b_factor: float
    element: str
    line: int

    @property
    def residue_id(self):
        return (self.chain, self.sequence, self.insertion, self.residue)

    @property
    def atom_id(self):
        return (*self.residue_id, self.name)

    def describe(self):
        return {'chain': self.chain, 'residue_number': self.sequence,
                'insertion_code': self.insertion, 'residue_name': self.residue,
                'atom_name': self.name, 'element': self.element, 'altloc': self.altloc,
                'occupancy': self.occupancy, 'b_factor': self.b_factor}


def parse_pdb(data: str):
    atoms = []
    issues = []
    current_model = 0
    model_records = False
    for number, line in enumerate(data.splitlines(), 1):
        record = line[:6].strip()
        if record == 'MODEL':
            model_records = True
            current_model += 1
            continue
        if record == 'ENDMDL':
            continue
        if record not in {'ATOM', 'HETATM'} or (model_records and current_model != 1):
            continue
        if len(line) < 54:
            issues.append(f'line {number}: coordinate record too short')
            continue
        try:
            def number_field(start, end, default):
                value = line[start:end].strip()
                if not value and default is None:
                    raise ValueError('missing coordinate')
                result = float(value) if value else default
                if not math.isfinite(result):
                    raise ValueError('non-finite value')
                return result
            element = line[76:78].strip().upper() if len(line) >= 78 else ''
            atoms.append(Atom(record, int(line[6:11]), line[12:16].strip(),
                line[16:17].strip(), line[17:20].strip().upper(), line[21:22].strip(),
                line[22:26].strip(), line[26:27].strip(),
                number_field(30, 38, None), number_field(38, 46, None), number_field(46, 54, None),
                number_field(54, 60, 1), number_field(60, 66, 0), element, number))
        except ValueError:
            issues.append(f'line {number}: invalid numeric field')
    # One conformer per atom; highest occupancy wins, with blank and A as tie breakers.
    selected = {}
    for atom in atoms:
        rank = (atom.occupancy, atom.altloc == '', atom.altloc == 'A')
        if atom.atom_id not in selected or rank > selected[atom.atom_id][0]:
            selected[atom.atom_id] = (rank, atom)
    return [item[1] for item in selected.values()], issues, current_model if model_records else 1


def analyze(data: bytes, ligand_name: str, chain: str | None, residue_number: str | None,
            insertion_code: str | None, cutoff: float, overlap: float):
    if not (math.isfinite(cutoff) and 0 < cutoff <= 20):
        raise ValueError('cutoff must be greater than 0 and at most 20 Å')
    if not (math.isfinite(overlap) and 0 <= overlap <= 3):
        raise ValueError('overlap must be between 0 and 3 Å')
    try:
        text = data.decode('utf-8')
    except UnicodeDecodeError as exc:
        raise ValueError('input is not UTF-8 text') from exc
    atoms, issues, model_count = parse_pdb(text)
    if not atoms:
        raise ValueError('no valid PDB coordinate records found')
    candidates = [a for a in atoms if a.record == 'HETATM' and a.residue == ligand_name.upper()
                  and (chain is None or a.chain == chain)
                  and (residue_number is None or a.sequence == residue_number)
                  and (insertion_code is None or a.insertion == insertion_code)]
    residues = sorted({a.residue_id for a in candidates})
    if not residues:
        raise ValueError('no matching HETATM ligand found')
    if len(residues) > 1:
        labels = ', '.join(f'{r[3]}:{r[0] or "(blank)"}:{r[1]}:{r[2] or "(blank)"}' for r in residues)
        raise ValueError('multiple ligand residues match; specify --chain, --residue-number and/or --insertion-code: ' + labels)
    ligand = [a for a in candidates if a.element not in {'H', 'D'}]
    receptor = [a for a in atoms if a.record == 'ATOM' and a.element not in {'H', 'D'}]
    if not ligand:
        raise ValueError('ligand has no non-hydrogen atoms')
    missing_elements = sum(1 for a in ligand + receptor if not a.element)
    unsupported_elements = sorted({a.element for a in ligand + receptor if a.element and a.element not in RADII})
    contacts = []
    for a in ligand:
        for b in receptor:
            distance = math.dist((a.x,a.y,a.z), (b.x,b.y,b.z))
            if distance <= cutoff:
                radii = RADII.get(a.element), RADII.get(b.element)
                clash = None if None in radii else (sum(radii) - distance >= overlap)
                contacts.append({'ligand_atom': a.describe(), 'receptor_atom': b.describe(),
                    'distance_angstrom': round(distance, 3), 'overlap_flag': clash})
    contacts.sort(key=lambda c: (c['distance_angstrom'], c['ligand_atom']['atom_name'],
                                 c['receptor_atom']['chain'], c['receptor_atom']['residue_number'],
                                 c['receptor_atom']['atom_name']))
    return {'schema_version': '1.0', 'tool': {'name': 'structural-audit', 'version': __version__},
        'input': {'format': 'PDB', 'sha256': hashlib.sha256(data).hexdigest(),
                  'models_seen': model_count, 'model_used': 1, 'parse_issues': issues},
        'selection': {'ligand': {'chain': residues[0][0], 'residue_number': residues[0][1],
                       'insertion_code': residues[0][2], 'residue_name': residues[0][3]},
                      'alternate_locations': 'highest occupancy per atom; blank then A on ties'},
        'parameters': {'contact_cutoff_angstrom': cutoff, 'overlap_threshold_angstrom': overlap},
        'summary': {'ligand_heavy_atom_count': len(ligand), 'receptor_heavy_atom_count': len(receptor),
                    'contact_count': len(contacts), 'overlap_flag_count': sum(c['overlap_flag'] is True for c in contacts),
                    'unclassified_contact_count': sum(c['overlap_flag'] is None for c in contacts),
                    'atoms_missing_element': missing_elements, 'unsupported_elements': unsupported_elements},
        'contacts': contacts,
        'limitations': ['First model only. ATOM records are treated as receptor; HETATM atoms are excluded from receptor contacts.',
                        'Distance-only overlap flags are triage signals; bonding, coordination and chemical context are not assessed.',
                        'Electron density, force fields, molecular dynamics and binding affinity are not evaluated.']}


def main(argv=None):
    p = argparse.ArgumentParser(description='Triage ligand–receptor contacts in a PDB structure.')
    p.add_argument('structure', type=Path, help='input PDB file')
    p.add_argument('--ligand', required=True, help='three-character PDB HETATM residue name')
    p.add_argument('--chain', help='ligand chain ID; use an empty string for a blank chain')
    p.add_argument('--residue-number', help='ligand residue sequence number')
    p.add_argument('--insertion-code', help='ligand insertion code; use an empty string for blank')
    p.add_argument('--cutoff', type=float, default=4.0, help='contact cutoff in Å (default: 4.0)')
    p.add_argument('--overlap', type=float, default=0.4, help='van der Waals overlap threshold in Å (default: 0.4)')
    p.add_argument('--output', type=Path, help='write JSON to this path; otherwise print to stdout')
    args = p.parse_args(argv)
    try:
        data = args.structure.read_bytes()
        report = analyze(data, args.ligand, args.chain, args.residue_number, args.insertion_code,
                         args.cutoff, args.overlap)
        result = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + '\n'
        if args.output:
            if args.output.resolve() == args.structure.resolve():
                raise ValueError('output must not overwrite the input structure')
            args.output.write_text(result, encoding='utf-8')
        else:
            print(result, end='')
        return 0
    except (OSError, ValueError) as exc:
        print(f'error: {exc}', file=sys.stderr)
        return 2
