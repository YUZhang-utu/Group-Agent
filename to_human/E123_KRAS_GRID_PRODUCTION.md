# Existing KRAS pools to docking MOL2

Run from the existing Group-Agent checkout used for E122. This add-on does not
modify the sealed prepared targets, original search code, E121 plan or exports.
No new whole-library coordinate search is performed. Compact catalog identity
arrays are scanned once; only pool members receive geometry evaluation.

Pull the production update in the same checkout used for the completed pilot:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin feature/structure-guided-chat
conda activate /mnt/medchem_taltio/wrk/yu_agent/Group_Agent/envs/aidd_workstation
export PYTHON="$(command -v python)"
nohup bash scripts/run_kras_grid_production.sh \
  "$HOME/kras-retained=v1/plan" "$HOME/kras-grid-production-v1" \
  > "$HOME/kras-grid-production-v1.log" 2>&1 < /dev/null &
tail -n 10 "$HOME/kras-grid-production-v1.log"
```

The equals sign in kras-retained=v1 is intentional. Twenty workers and one
numeric-library thread per worker are defaults. The launcher uses Linux flock
to prevent simultaneous runs on the same output. Use the wrapper for execution.
No tmux is required. Keep the original library/source files available throughout.

Stages: verify source seals and native controls; map at most three conformers per
pool ID; score with the unchanged E122 grid; write tier and failure manifests;
source-check and export G12D; source-check and export G12V. Source IDs are mapped
using compact arrays, not by redoing the original full-library pose search.
Original raw MOL2 files are streamed to retrieve selected records during export.
Pilot results estimate approximately nine hours for pool scoring alone, excluding
mapping, checks and export; this is not a same-day completion guarantee.

Selection: balanced-profile passes, ordered by that profile's supporting pose
score, up to 1,000,000 unique library molecule IDs per target. No random padding,
activity claim or chemical-diversity claim. All strict passes are inside balanced,
but strict is not a separate score-order preference. With pilot rates, expect
approximately 255,000 G12D candidates and a budgeted 1,000,000 G12V candidates
before source checks; actual counts may differ. No replacement is silently added
for quarantined source records. This is an explicit docking budget, not proof
that non-exported molecules are unsuitable.

Results:

- report.json: running stage or complete status, exact counts and export locations.
- results.sqlite: committed results for every original target/molecule member,
  with profile flags and supporting transforms. No pool IDs are discarded.
- all_molecule_tiers.csv: balanced, lenient_reserve or unresolved for every ID.
  Balanced includes IDs beyond the export quota; molecules.csv identifies exports.
- failure_diagnostics.json: overlapping failures of the best exact scored pose
  for non-balanced molecules; not an exhaustive analysis of all poses.
- KRAS_G12D_export.json and KRAS_G12V_export.json: completed target receipts.
- exports/KRAS_G12D_attempt_001/mol2/batch_001.mol2 (and analogous G12V):
  at most 100,000 molecules per MOL2, final file may be smaller.
- Each attempt's molecules.csv maps IDs/conformers to batch and record;
  quarantined.csv lists source failures. Quarantine is not biological rejection.

The selected conformer's original coordinates, bond records and atom information
are preserved. These are docking inputs, not receptor-aligned predicted poses.
Normal docking preparation (protonation, charges and format requirements of your
docking software) is still needed. Even profile-passing poses can have clashes.

To resume after a stopped process, run the same command with the same paths and
settings. Committed scoring batches are skipped. A completed target export is
hash-verified and reused. An interrupted target export starts a new numbered
attempt; old partial files stay on disk. Use ONLY the folder identified by its
completed target receipt, not an older partial attempt. Keep sufficient disk
space for the database, WAL, MOL2 files and any interrupted export attempt.
The original E121 exports retain the reserve/unresolved structures for later work.

Validation: nine focused local tests passed, including native controls, pool-only
mapping, resume, input-change rejection, source-coordinate mismatch quarantine,
native source geometry checking and MOL2 splitting. Test geometry and I/O fixtures
are separate; these tests do not validate activity or full-scale throughput.
English-content and shell syntax checks passed. School-workstation execution is
pending; no actual production molecules have been generated on that host yet.
