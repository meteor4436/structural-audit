# Structural Audit

**A structure is a claim. Make the evidence inspectable.**

Structural Audit is a small, reproducible starting point for reviewing a protein–ligand model. It reads a PDB coordinate file, selects one ligand residue, enumerates nearby receptor atoms and flags unusually close contacts. Every result carries the input SHA-256, atom identities, selection rules and thresholds in machine-readable JSON.

This repository implements **coordinate-based triage**. It is deliberately narrow: a close contact is a prompt for chemical and structural review, not evidence of misconduct or proof that a ligand is wrong.

## Run in 30 seconds

Python 3.10+ is the only runtime requirement. From the repository root:

```bash
python3 -m structural_audit examples/synthetic.pdb --ligand LIG --chain B --residue-number 99
```

The example is **synthetic**, designed to exercise the checker rather than represent an experimental structure. Its output contains one carbon–carbon contact at 2.7 Å and one overlap flag.

For your own PDB file:

```bash
python3 -m structural_audit structure.pdb --ligand LIG --chain A --residue-number 501 --output audit.json
```

If a ligand name matches more than one residue, the command stops and asks for a more specific selection. Use `--insertion-code` when residue numbering alone is ambiguous. Blank chain or insertion code can be specified with an empty string, e.g. `--chain ''`.

Optionally install the command line entry point:

```bash
python3 -m pip install .
structural-audit examples/synthetic.pdb --ligand LIG --chain B --residue-number 99
```

## What the report contains

- The input file's SHA-256 and the model used, to identify exactly what was analyzed.
- Ligand chain, residue number, insertion code and atom identity.
- Sorted atom-pair distances up to a configurable cutoff (default 4.0 Å).
- A distance-only van der Waals overlap flag (default overlap ≥ 0.4 Å), using a small explicit element-radius table.
- Parse issues, unsupported or missing element fields, unclassified contacts and analysis limitations.

Only the first PDB model is used. Alternative locations are selected per atom by highest occupancy, preferring blank and then `A` on ties. This rule is reproducible but can combine alternate conformations; inspect the source structure before interpreting the result. Only `ATOM` records are treated as receptor atoms; water, ions, cofactors and other `HETATM` records are excluded. Hydrogens and deuterium are excluded. Unrecognized elements do not receive an overlap classification. Covalent ligands, metal coordination and bonded topology require human review.

```bash
python3 -m structural_audit --help
python3 -m unittest discover -s tests -v
```

## Research direction

Future work could add mmCIF input, model-to-map fitting and independent half-map checks, chemically aware contacts, parameterized molecular dynamics, and reproducible free-energy calculations. None of these analyses are implemented or validated here. Any future claim about density, affinity or model quality should include the original data, method, uncertainty and independent review.

## License

MIT. See [LICENSE](LICENSE).
