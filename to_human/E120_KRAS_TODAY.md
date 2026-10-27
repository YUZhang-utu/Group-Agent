# E120: recover physically plausible contact poses from the existing direct pool

This replaces the zero-yield E118 launch, without changing any mandatory contact.
The input molecule IDs, source conformers, reference packages and sealed original
code remain unchanged. No full-library counting job is repeated.

## Run

Activate the workstation environment, pull `feature/structure-guided-chat`, then:

```bash
export PYTHON="$(command -v python)"
nohup bash scripts/run_kras_refine_today.sh \
  "$HOME/kras-direct-count-v1" \
  "$HOME/kras-soft-clash-v1/report.json" \
  "$HOME/kras-refine-v1" 4 \
  > "$HOME/kras-refine-v1.log" 2>&1 < /dev/null &
```

The final `4` is a cooperative batch compute budget in hours. Startup validation,
the initial saved-pose gate, draining in-flight chunks and source export are extra.
It is not a guaranteed four-hour total runtime or a million-molecule yield promise.
Use a smaller positive number if docking must start sooner.

## Automatic stages

1. Read the E119 saved examples (up to 20 per target), verify molecule/conformer
   identity and replay their original coordinates, distance and required contacts.
2. Apply rigid optimization with simultaneous receptor-repulsion and all required
   contact penalties. No internal molecular distortion, torsion edits or receptor
   relaxation. Independent final acceptance retains all original contact masks
   and minimum receptor-heavy-atom distance >=2.0 A. Perturbed native controls
   must also pass. If none of the actual saved poses can be rescued, stop rather
   than launch an expensive unproductive batch.
3. Warm-start a resumable result database from rescued examples. Process a spread
   order of existing conformers, skipping IDs outside the existing direct pool.
   Generate finite initial poses only for pool molecules, then refine up to three
   low-overlap contact-passing starts. This recovers poses not saved by the old
   counter. Failed optimization does not establish that a molecule cannot bind.
4. Stop at the compute budget, plan completion, or the requested per-target quota;
   drain active chunks and preserve all committed results. Report measured yield
   and coverage. A partial run is explicitly labelled as budgeted, not full coverage.
5. Select up to 1000000 eligible molecule IDs per target using the prior 80/20
   score/diversity-proxy policy and export at most 10 MOL2 files per target, up to
   100000 molecules per file. A final smaller file and shortfall are explicit.
   Never add clashing molecules to fill a quota. The export path is timestamped.

Export verifies each original MOL2 record hash, parsing, atom order and source
coordinates against the artifact, then rechecks contacts and minimum distance
using exact original coordinates under the saved pose transform. Failures are
listed in `rejected.csv`. MOL2 output preserves source coordinates/chemistry and
changes only the title. The evaluated transform is in the molecule CSV. One best
scored source conformer is exported per library molecule ID; canonical chemical
deduplication and flexible macrocycle docking preparation are separate.

## Status and files

```bash
tail -n 30 "$HOME/kras-refine-v1.log"
cat "$HOME/kras-refine-v1/gate/report.json"
cat "$HOME/kras-refine-v1/batch/progress.json"
```

The gate report contains real selected-example rescue counts, not an unbiased
full-pool yield estimate. Batch progress reports actual eligible molecule counts.
The final log prints `export-TIMESTAMP/report.json`. Each target folder contains
`batch_001.mol2` onward and `molecules.csv`. Check `quota_met`, target `shortfall`,
source coverage and `exact_source_pose_rechecked` before handoff.

To add another compute budget with the same code and inputs, run the same command
again. Committed chunks are skipped, and a fresh export snapshot is created.
Do not run two instances against the same batch directory or update code/input
files mid-run. A failed export does not discard refinement results; rerun the
export CLI against `batch/` with `--allow-budgeted` and a new output path.

## Limits

This is contact-constrained rigid pose recovery and exploratory ranking, not
validated affinity prediction or finished docking. The retained distance floor
is a gross-clash proxy, not a complete atom-specific VDW/energy test. Required
polar proximity is not proof of proper hydrogen-bond angles. Shape, directional
feature matching and reference-region coverage remain soft scores; whole-ligand
color is diagnostic only. Conformer sampling, rigid bounds and the time budget
limit recall. The actual candidate yield is unknown until the workstation runs.
