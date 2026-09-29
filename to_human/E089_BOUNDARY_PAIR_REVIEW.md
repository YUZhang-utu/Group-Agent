# Boundary versus definite-class pair review

This stage computes actual source-backed geometry and continuous side-chain
property/steric comparisons. It does not change memberships or adopt a cutoff.
It runs using the completed backbone build and E088 candidate database, without
requiring complete chemistry/typed builds.

## Workstation command

Update feature/structure-guided-chat and run from the repository root:

```bash
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
BUILD=/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/backbone
python -m aidd_agent.boundary_pair_review \
  --build "$BUILD" \
  --candidates "${BUILD}-boundary-candidates-v1" \
  --output "${BUILD}-boundary-pair-review-v1" \
  --per-band 100 --references 3
```

Use a fresh output. Run in the RDKit environment used for the completed build.
The command verifies the candidate database hash and input receipts. It samples
up to 100 original classes in each angular band: <=35, 35-45 and >45 control.
Each sampled class contributes one minimum-hash conformer per band. All proposed
target classes/rotations are retained. Each target supplies up to three distinct
molecule references (first indexed records, not population-unbiased samples).

Only selected source records undergo chemistry reconstruction. Needed MOL2 files
are each scanned once until all selected records are found. This is not a full
library chemistry rerun, but sparse access can still require substantial text I/O.
Progress reports include target selection, source files and record locations.

For every selected record, the raw hash, names, original descriptor values and
stored residue atom mapping must reproduce. A changed or missing record stops
the run without a complete report. New typed property ordering is explicitly
aligned to the stored backbone ordering before candidate rotations are applied.

## Read these outputs

- report.json: sampled counts, missing within-class controls and metric summaries.
- pair_metrics.csv: every candidate/reference and within-definite reference pair.
- panel.json: exact selected IDs, target classes and directed rotations.
- profiles.jsonl: per-record provenance, continuous properties and steric proxies.
- selected.mol2: original selected coordinates, not receptor-prepared docking input.

Metrics remain separate:

1. Backbone RMSD after proper rigid alignment (angstroms; no reflections).
2. Scaled continuous side-chain property RMS difference.
3. Scaled local steric-profile RMS difference.
4. Maximum per-residue steric-profile difference, to expose a local outlier that
   could be hidden by an average.

Compare candidates against controls from the same target class. Do not use a
single pooled median as a merge cutoff: targets contribute unequal numbers of
pairs, and a few first-indexed references do not describe all of a large class.
Classes with only one distinct reference molecule have no within-class control
and are explicitly counted. No minimum-score cherry-picking is performed.

This panel evaluates plausibility, not final full-library acceptance or search
recall. It does not automatically export a complete sample of every class,
reassign any boundary record, run PLANTS or run N-E rescoring.

## Local validation

Seven focused tests passed (selection, proper-rotation RMSD, explicit unit
mapping, candidate rotation direction and prior angle/property regressions).
Existing small dataset completed source-backed extraction: 3 boundary records,
3 reference records, 5 boundary/reference pairs and 1 within-class control.
One target lacked enough distinct molecules for a control. This is fixture
evidence only; workstation measurements remain pending.
