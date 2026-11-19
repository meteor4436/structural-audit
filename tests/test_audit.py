import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from structural_audit.cli import analyze, parse_pdb


def pdb_line(record, serial, name, residue, chain, seq, x, y, z, element, alt='', occupancy=1.0, insertion=''):
    return (f'{record:<6}{serial:>5} {name:>4}{alt:1}{residue:>3} {chain:1}{seq:>4}{insertion:1}   '
            f'{x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{20.0:>6.2f}'
            f'          {element:>2}  \n')


class StructuralAuditTest(unittest.TestCase):
    def fixture(self):
        return (pdb_line('ATOM', 1, 'CA', 'ALA', 'A', 7, 0, 0, 0, 'C') +
                pdb_line('ATOM', 2, 'N', 'ALA', 'A', 7, 12, 0, 0, 'N') +
                pdb_line('HETATM', 3, 'C1', 'LIG', 'B', 99, 2.7, 0, 0, 'C') +
                pdb_line('HETATM', 4, 'O', 'HOH', 'A', 8, 2.8, 0, 0, 'O')).encode()

    def test_contact_clash_and_provenance(self):
        report = analyze(self.fixture(), 'LIG', 'B', '99', None, 4, 0.4)
        self.assertEqual(report['summary']['contact_count'], 1)
        self.assertEqual(report['summary']['overlap_flag_count'], 1)
        self.assertEqual(report['contacts'][0]['distance_angstrom'], 2.7)
        self.assertEqual(report['contacts'][0]['receptor_atom']['residue_name'], 'ALA')
        self.assertEqual(len(report['input']['sha256']), 64)

    def test_multiple_ligands_need_disambiguation(self):
        data = self.fixture() + pdb_line('HETATM', 5, 'C1', 'LIG', 'C', 100, 1, 0, 0, 'C').encode()
        with self.assertRaisesRegex(ValueError, 'multiple ligand residues'):
            analyze(data, 'LIG', None, None, None, 4, 0.4)
        self.assertEqual(analyze(data, 'LIG', 'C', None, None, 4, 0.4)['selection']['ligand']['chain'], 'C')

    def test_altloc_uses_highest_occupancy_and_first_model(self):
        data = ('MODEL        1\n' + pdb_line('HETATM', 1, 'C1', 'LIG', 'A', 1, 1, 0, 0, 'C', 'A', 0.3) +
                pdb_line('HETATM', 2, 'C1', 'LIG', 'A', 1, 3, 0, 0, 'C', 'B', 0.7) +
                'ENDMDL\nMODEL        2\n' + pdb_line('HETATM', 3, 'C1', 'LIG', 'A', 1, 9, 0, 0, 'C') + 'ENDMDL\n')
        atoms, errors, count = parse_pdb(data)
        self.assertEqual((len(atoms), atoms[0].x, count, errors), (1, 3, 2, []))

    def test_unknown_element_is_not_classified_as_clash(self):
        data = self.fixture() + pdb_line('ATOM', 8, 'ZN', 'HIS', 'A', 8, 3, 0, 0, 'ZN').encode()
        report = analyze(data, 'LIG', 'B', '99', None, 4, 0.4)
        self.assertEqual(report['summary']['unclassified_contact_count'], 1)
        self.assertEqual(report['summary']['unsupported_elements'], ['ZN'])

    def test_cli_json_and_errors(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'sample.pdb'
            path.write_bytes(self.fixture())
            run = subprocess.run([sys.executable, '-m', 'structural_audit', str(path), '--ligand', 'LIG'],
                                 text=True, capture_output=True, check=True)
            self.assertEqual(json.loads(run.stdout)['summary']['contact_count'], 1)
            fail = subprocess.run([sys.executable, '-m', 'structural_audit', str(path), '--ligand', 'BAD'],
                                  text=True, capture_output=True)
            self.assertEqual(fail.returncode, 2)
            self.assertIn('no matching', fail.stderr)
            same = subprocess.run([sys.executable, '-m', 'structural_audit', str(path), '--ligand', 'LIG', '--output', str(path)],
                                  text=True, capture_output=True)
            self.assertEqual(same.returncode, 2)
            self.assertEqual(path.read_bytes(), self.fixture())

    def test_missing_coordinates_reported(self):
        bad = pdb_line('ATOM', 1, 'CA', 'ALA', 'A', 7, 0, 0, 0, 'C')
        bad = bad[:30] + '        ' + bad[38:]
        atoms, errors, _ = parse_pdb(bad)
        self.assertEqual(atoms, [])
        self.assertEqual(len(errors), 1)


if __name__ == '__main__':
    unittest.main()
