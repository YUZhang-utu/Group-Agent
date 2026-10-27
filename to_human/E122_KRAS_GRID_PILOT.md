# Fast receptor-grid pilot: 7RPZ / 12FF

Use the E121 plan, not the export folder. The returned report shows the actual
workstation directory is `kras-retained=v1` (equals sign), not `kras-retained-v1`.
Keep that exact path below unless your actual plan is elsewhere.

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
conda activate /mnt/medchem_taltio/wrk/yu_agent/Group_Agent/envs/aidd_workstation
export PYTHON="$(command -v python)"
git pull --ff-only origin feature/structure-guided-chat
nohup bash scripts/run_kras_grid_pilot.sh \
  "$HOME/kras-retained=v1/plan" "$HOME/kras-grid-pilot-v1" 10000 \
  > "$HOME/kras-grid-pilot-v1.log" 2>&1 < /dev/null &
tail -n 10 "$HOME/kras-grid-pilot-v1.log"
```

Stage sequence: seal checks and receptor grid controls; sample ID mapping;
grid_scoring; report. Mapping reads only compact catalog identity arrays, once;
scientific evaluation loads only sampled pool conformers. No full-library pose
search, torsion optimization, MOL2 export or automatic follow-on job occurs.
Use a fresh output directory for each run. Pilot output is not resumable.

The sample is 10,000 unique pool molecule IDs per target, shared when possible,
with up to three existing conformers per ID. Conformer selection is deterministic
by global ID, not an assertion that these are the most relevant conformers.
Original finite 517-ish seed generation and mandatory contact rules remain.
Grid scores shortlist eight poses per conformer for exact distance/burial checks.
Out-of-grid poses remain unresolved. Reference-ligand coordinates locate the
box only; burial is computed from protein occlusion, not reference shape overlap.
GDP/Mg are included as in the sealed receptor preparation; waters are absent.

Read `$HOME/kras-grid-pilot-v1/report.json` when status is complete. It reports:

- Native control retention and interpolation error in controls.json.
- Each target's sampled molecules, actual conformers, missing mappings and direct
  contact recovery. No-recovered-contact molecules remain unresolved.
- Lenient, balanced and strict exploratory profile retention, projected counts
  and sampling intervals. Profiles allow some overlaps; they are not valid-pose
  certification. A failed profile does not prove molecular impossibility.
- Compute time separately from setup/index time; rough full-pool compute estimate
  excludes export and may change with sharing and library composition.
- Per-molecule outcomes and supporting transforms in molecules.jsonl; best score
  and each passing profile's own supporting pose are recorded separately.

All profiles must be interpreted with their native-control result. Two native
self-controls establish basic consistency, not enrichment or false-negative rate.
Do not silently tune profiles until a million-ID quota is met. If more than a
million remain, later ranking can implement an explicit work budget while retaining
the rest as unresolved. This pilot creates no production selection and leaves all
original pools and E121 exports unchanged.
