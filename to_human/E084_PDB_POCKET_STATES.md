# E084: experimental PDB pocket states

Independent export review found a critical geometry issue in the original C5
representative, 7BJ0:A: the fitted core RMSD is 1.82 A, but all observed pocket CA
atoms have RMSD 23.77 A and maximum displacement 45.80 A. Do not treat this state
as an accepted consensus/docking receptor. The original six-cluster result is
preserved for analysis; the all-pocket alignment guard requires follow-up before
production acceptance. See the E084_MDM2_ANALYSIS_REVIEW package, tables 10 and 11.

The first implementation uses the existing verified diversity collection as input,
but compares receptor pocket space rather than ligand fingerprints. New diversity
reports require this review before consensus; old sealed reports remain legacy.

## Computation and permissive defaults

- Align canonical-residue CA coordinates using a recorded, trimmed local core.
- Build one common 1-Angstrom grid around the reference ligand, expanded by 4 A.
  Remove receptor VDW spheres plus a 0.5-A probe. This is a bounded local cavity
  descriptor, not full cavity discovery or a ligand entrance/kinetic simulation.
- Compare cavity IoU and local subregions. Geometry contributes 85%, chemistry 15%.
  Local differences require at least 30 cubic Angstroms of changed space before
  contributing a supplemental local term; tiny variations do not each define states.
- Six residue-template channels: donor, acceptor, hydrophobic, aromatic, positive,
  negative. Donor/acceptor positions use an outward heavy-atom direction proxy.
  These annotations are not protonation calculations, electrostatic energies,
  hydrogen reconstruction or experimentally validated interaction strengths.
- Weighted average-linkage clustering uses distance 0.40 by default; 0.30 and 0.50 are
  reported for sensitivity only. This composite distance is not an IoU cutoff.
  There is no forced cluster count and no automatic small-cluster deletion.
- Select a real quality-aware medoid. Each PDB contributes total unit weight within
  a cluster regardless of chain multiplicity. Linkage also shares each PDB's unit
  weight across its candidate chains. PDB support is not state occupancy.
- Flag broad clusters whose medoid does not cover every member within the cutoff,
  and substantial within-cluster chemical differences. Inspect before adoption.

Missing/ambiguous atoms near the reference pocket cannot be treated as extra
space. Missing remote sidechain atoms outside the reference exclusion envelope
do not alone disqualify a pocket. The atom scope is recorded. This reference-based
rule cannot rule out an unobserved remote sidechain moving into the pocket.
Identity assembly 1, mapped target chains and model 1 are the supported scope.
Other chains/interfaces, metals, waters, nonstandard residues and unresolved
occupancy/protonation require review. Pocket mutations are held separately.

## Chat sequence

1. Collect structures with the verified structure_diversity workflow.
2. Request pocket analysis (`pockets` intent), specifying actual reference_query
   and target_chain if ambiguous. Changing cluster_distance creates a fresh run.
3. Inspect report.html, overlap/distance matrices, sensitivity and held structures.
4. Explicitly adopt all states, or adopt with pocket_selection.cluster_ids for a
   subset. Automated chains stop here; general continuation is not approval.
5. Request consensus with reference.pocket_state_id from the adopted report.
   One state per task; no incompatible contacts are pooled across states.
6. Continue existing recommendation/design/search workflows on that state.

Natural language is routed to these intents. Do not assume every front-end has
dedicated slash-command parsing for the new intents. A target alone does not
authorize the agent to invent a reference binding site when multiple sites exist.

## Direct commands

With the existing AIDD environment and repository src on PYTHONPATH:

```bash
python -m aidd_agent.pocket_states build --source /path/to/diversity/report.json --output /new/pockets --reference-query 5C5A:NUT:A:201 --target-chain A
```

The example reference is MDM2-specific. Use verified target-specific IDs elsewhere.
After user review, explicitly adopt:

```bash
python -m aidd_agent.pocket_states adopt --source /new/pockets/report.json --output /new/pockets-adopted
python -m aidd_agent.pocket_states consensus --source /new/pockets-adopted/report.json --output /new/state-consensus --pocket-state-id ACTUAL_ADOPTED_ID
```

The consensus stage uses existing CCD preparation and may fetch ligand CCD data.
The adoption step preserves source/report hashes. Changed structures, grids or
analysis reports invalidate consensus dispatch. Rebuild and review changed inputs.

## Evidence and limitations

The local MDM2 experiment uses the existing 133-file E052 collection and one
reference region, not a new exhaustive RCSB census. Inspect the E084 validation
receipt for actual accepted/held chains and cluster counts. Initial settings
produced excessive subdivision; the exploratory default was relaxed and local
quality scope corrected, with sensitivity retained. This is development-set
calibration, not an independent docking-recall benchmark.

Representative structures have not yet been shown to preserve every ligand pose
or scaffold. In particular, membership in a broad cluster does not prove that
every local MDM2 opening is interchangeable with its medoid. No PLANTS or N-E run
is claimed by this feature. Those remain later validation stages.
