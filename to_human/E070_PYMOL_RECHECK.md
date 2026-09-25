# E070: interaction recheck on real MDM2 coordinates

The earlier tests validated synthetic detection and mocked drawing, not chemical
coverage on the user's complexes. Two functional gaps were identified: only
ChemPy bond order 4 contributed aromatic atoms, and hydrophobic contacts were
excluded by default. Deposited mmCIF aromatic flags now supplement live bond
typing; deduplicated hydrophobic contacts are included by default. Nearby protein
selection retains whole residues to avoid truncating ring geometry at the cutoff.

Exploratory offline replay, 2026-09-25, unchanged geometry thresholds:

| Complex | Pi candidates | Hydrophobic atom pairs | Hydrophobic residue representatives |
|---|---:|---:|---:|
| 6Q9L:HTZ:A:201 | 1 | 12 | 5 |
| 3JZK:YIN:A:1 | 1 | 11 | 5 |

Raw replay reports are in `data/e070-6Q9L/report.json` and
`data/e070-3JZK/report.json` on the development machine. Their input SHA256 values:

- 6Q9L: `71071f1999e716797762e05a0a90981a5b4d2a5cca968edac864921c7928bf61`
- 3JZK: `8188cda0b95ef009cf2f56a3e043e458cdea480e7a077e051cf0872069fd456c`

This replay uses actual crystal coordinates and deposited component bonds. It is
not a PyMOL rendering test, does not evaluate polar/halogen acceptor typing, and
does not prove that the user's live scene has the same bonds or coordinates.
No compatible PyMOL package was available in the development Windows runtime.
Salt/cation-pi remain limited to loaded formal charges, with visible coverage
warnings; protonation, water bridges and metal coordination remain unimplemented.

## Workstation verification

Pull the feature branch, restart Chat and its bridge PyMOL, then reopen the
completed consensus task with `/pymol TASK_ID`. Do not rerun the scientific task.
Use this direct request first to bypass LLM routing:

```text
/view {"operation":"typed_interactions","types":["polar_contact","pi_stacking","cation_pi","salt_bridge","halogen_bond","hydrophobic"]}
```

Inspect `/pymol_status`: the receipt now includes `chemistry`,
`chemistry_warnings`, and `type_counts` (accepted/displayed/requested for all six
types). Objects are named `vNNN_ix_pi_stacking` (cyan),
`vNNN_ix_hydrophobic` (gray70), and so on. A zero count must not be replaced by
invented lines. Check `component_status=available` and ring counts when pi
assignments are missing. The JSON/CSV links retain individual atom/residue evidence.

For reproducible real-PyMOL validation, run in a runtime with PyMOL, NumPy and Gemmi,
from the repository root. Set `STRUCTURE` to the existing 6Q9L mmCIF path:

```bash
PYTHONPATH=src python scripts/check_pymol_interactions.py \
  --structure "$STRUCTURE" --ligand HTZ --chain A --resi 201 \
  --output pymol-check-6Q9L
```

The output directory must be new. This runs detection twice, checks measurement
objects and dash colors, checks repeatability, and writes PNG/PSE and reports.
Add `--offline` only for coordinate-only validation without PyMOL. Natural-language
verification follows the successful direct request using the prompt in E069.
