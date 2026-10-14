# Screening components and benchmark readiness, 2026-10-14

## Evidence boundary

Inspected current source, historical workflow records and returned E108 outputs in
D:/agent/aidd, including historical_recovery.csv. No live workstation query NPZ,
adopted-design/report.json or old/current protocol.json has been supplied locally.
The twelve returned query hashes agree, but do not reveal their settings or prove
equality with the old query. Local E053/E054/E056 preparation reports are historical
development artifacts and must not be substituted for the running adopted query.

## Component audit

| Component | Implementation / current execution evidence | Remaining check |
|---|---|---|
| Crystal census and quality | E052 census, target/site/occupancy/completeness checks | Exact source cohort of the adopted query |
| Ligand diversity | Morgan radius 2, 2048 bits, chirality; farthest-first chemical coverage | Exact eleven adopted IDs, coverage and exclusions; not 3D diversity by itself |
| Protein alignment/state admission | Local pocket correspondence, CA RMSD and severe-overlap admission | Actual cohort/reference and excluded alternative states |
| Contact evidence | Protein residue/atom-specific spatial/directional consensus | Selected anchor IDs, evidence, mandatory versus optional roles |
| Coarse filtering | Current threshold compute enforces adaptive_coarse when present; otherwise coarse_constraints | Actual limits; these are alternative branches, not necessarily two active sequential gates |
| Anchor geometry | Threshold path applies required/alternative constraints | Empty required lists do not create mandatory interactions; inspect adopted values |
| Shape/color matching | same_pose_gaussian uses 0.5 shape plus 0.5 color Tanimoto | Query template arrays, coordinate frame and Gaussian weight |
| Receptor exclusions | pocket_mask checks transformed ligand points against one receptor resource plus declared exclusions | Actual receptor identity, cutoff, clash fraction, exclusions and reference controls |
| Composite pose acceptance | Threshold path enforces minimum_pose_score | Numeric threshold, weights, calibration provenance |
| Contact-first output ranking | Current per-template contact/Gaussian order, positive-contact tier, RRF k=60 | Preserve this policy in both baseline and block experiment |
| Spatial occupancy rewards | Implemented optional spatial_groups/occupancy_weight | Not established as enabled; no automatic subpocket occupancy requirements |
| Pocket-state ensemble | E084/E093 state analysis/adoption supported for new projects | Legacy query not retrofitted; one receptor resource is not an ensemble |
| Water/metal-mediated features | Evidence retained but selected mediator anchors are rejected by contact_policy | Do not turn unsupported mediation into active screening constraints |
| Block priority | Twelve E106 arms, ChemPLP/EquiScore, receptor context 7NA2_A | Compatibility with search reference and whether priority predicts reference recovery |
| Stop rule | Complete a block, stop after 100000 qualified conformers or selected-block exhaustion | This does not certify global Top-K recovery |

Gaussian evaluation is used by the current engine; its composite contribution can
be zero if gaussian_weight is zero. Gaussian also has a separate contact-tie role
in final ranking. Presence of a module does not establish a nonzero runtime weight.
The shape reference is a ligand, not the receptor cavity. Pocket collision masks
are not full energetic complementarity, docking optimization or induced fit.

## Actual-configuration collection

scripts/audit_screening_configuration.py is standalone, standard-library only.
After pulling the delivered update on the workstation, run:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin feature/structure-guided-chat
python scripts/audit_screening_configuration.py --output "$HOME/mdm2-screening-config-audit"
```

Alternatively copy the standalone script and run it directly. No installation is
required. This command does not change or restart existing scientific jobs.

Defaults locate the original eleven-template search and all twelve E106 runs.
Use --runs-root and --historical-search if directories moved. Use a fresh output.
It reads sealed current query/protocol files, extracts scientific settings, hashes
declared receptor/template resources, compares old/current settings and records gaps.
It does not submit search/docking, alter parameters, or write to scientific inputs.
Old report/protocol files are fingerprinted but not externally authenticated by a
parent seal in this tool. Missing or unsealed resource provenance remains explicit.

Return the entire small JSON/report output folder to D:/agent/aidd. The audit includes
scientific protocol paths and code hashes, not environment variables or credentials.
Do not treat equal settings as equal execution: budget versus threshold branching
must be compared separately. Keep the original files for source verification.

## Benchmark readiness decision

Not ready to assert that every supported component is enabled or that all relevant
physics is modeled. The essential structural screening components already exist;
the immediate problem is runtime-policy consistency and validated routing, not
blindly enabling every module. Optional features must have compatible evidence.

Before performance experiments:

1. Recover the actual configuration and identify query reference versus 7NA2_A.
2. Freeze one adopted evaluator, template panel, receptor state, geometric gates,
   score weights, pose search, ranking, eligibility and output budget.
3. On molecule-disjoint fixed test candidates, compare candidate-level pass/fail,
   raw score and pose outputs with/without block scheduling. With all eligible
   blocks visited, global aggregation must reproduce the common full evaluator.
4. Only then measure priority-order curves and early stopping versus random blocks.
   Do not use a different evaluator or global-reference oracle for production stops.
5. Test nested library scales and multiple independent queries, reporting construction,
   target-specific priority scoring, online search, cache policy and memory separately.

Configuration/equivalence audits are acceptance checks, not biological validation.
Current low historical recall is dominated by selected/visited-block coverage, but
the historical ANN reference is not the same evaluator's exhaustive ground truth.

## Local validation

Four new synthetic configuration-audit tests plus fourteen E108 tests passed.
English-content guard passed. No actual workstation configuration was inspected
and no scientific execution or GitHub delivery was performed in this audit turn.
