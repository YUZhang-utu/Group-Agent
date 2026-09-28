# E082: MDM2 block panels, PLANTS receptor ensembles and N-E rescoring

Design protocol, 2026-09-28. No docking or inference has been run under this
protocol. User decisions: PLANTS, MDM2, existing in-house N-E workflow for
rescoring; binding-affinity models remain undecided and optional.

## Hypothesis

Combining representative block panels with receptor-state diversity and N-E
rescoring may recover useful molecules missed by one co-crystal/3D reference.
Test this against single-receptor and 3D-only baselines at equal measured budgets.
Predictions are not experimental affinity evidence.

## Keep sampling independent of the crystallographic ligand

Start with 100 unique molecules per validated conformer block, including their
actual in-block conformers. Expand promising or uncertain blocks to a nested 200;
small blocks use all molecules. Evaluate the SAME panel through existing 3D
search and docking. Do not require a good 3D score before docking the initial
representative panel. Otherwise that filter can irreversibly remove different
scaffolds before any ensemble benefit is possible.

The existing sampler supports this without changing production defaults:

```bash
python -m aidd_agent.conformer_block_store \
  --model /path/to/validated/typed/capacity-10000 \
  --output /path/to/new-panel-100-200 --sample-sizes 100 200 --seed 20260928
```

This generates membership manifests only, not dock-ready structures or engine jobs.

## MDM2 receptor-state panel

Reuse the existing Q00987-verified experimental structure census and aligned
MDM2 cohort. Select representative RECEPTOR pocket states, not merely diverse
co-crystal ligand fingerprints. Inspect local backbone geometry, Tyr100 and His96
rotamers, pocket shape/volume, lid-region coverage and construct differences.
Start with approximately 4-8 experimentally supported representatives if distinct
qualified states exist; this is a pilot budget proposal, not an established optimum.
Include peptide-bound states when their target identity and pocket preparation pass.
Do not count duplicate crystal copies as independent state support.

Check sequence, mutations, missing residues, alternate locations, protonation,
waters/cofactors and alignment. Distinguish an unresolved/truncated lid from a
measured open lid. Define comparable spatial search regions across aligned states
using the pocket ensemble, rather than a tight envelope of one small ligand.
Keep every receptor as a physically consistent structure; do not union all receptor
atoms into an artificial single pocket or demand all mutually exclusive contacts.

Freeze the SAME receptor panel and state count for every sampled ligand. PLANTS
ensemble docking approximates receptor flexibility by evaluating different
structures; it is not by itself continuous induced-fit simulation. Confirm the
installed PLANTS version/manual and actually supported ligand/ring/side-chain
flexibility before constructing executable settings. No PLANTS parameter syntax
or N-E interface is inferred in this protocol.

## Pose and block attribution

Retain source_block_id and source_conformer_id. If ligand ring geometry changes
during docking or N-E refinement, also retain output pose geometry, parent pose,
and a recomputed pose_block_id or out_of_domain flag. A favorable migrated pose
does not establish that the original frozen conformer geometry was favorable.
Preserve chemical identity and separately track protonation/microstates.

Rescore each pose with the SAME receptor state used for its generation, unless
the N-E workflow explicitly produces a new receptor/pose pair with new hashes.
Do not mix the best docking pose from state A with an unrelated N-E score from
state B and call them a single jointly supported pose.

## N-E interface and optional affinity

Use the user's N-E workflow; do not expand or invent the meaning of its name.
Its callable entry point, input/output example, score direction, units, model
version and whether it changes coordinates must be supplied or located before
an executable adapter can be validated. Future affinity plugins are optional;
their absence is not failure and their outputs must not be fabricated.

The result table preserves at least:

`run_id, target_id, receptor_state_id, receptor_sha256, protocol_id, signal_family,
model_name, model_version, molecule_id, source_conformer_id, source_block_id,
microstate_id, pose_id, parent_pose_id, pose_sha256, pose_block_id, status,
raw_score, score_unit, better_direction, pose_qc_pass, elapsed_seconds`.

Statuses distinguish ok, preparation_failed, docking_failed, no_valid_pose,
unsupported_input, model_failed and not_run. Missing results are not zero affinity.
If a future model only takes molecular/sequence inputs, cache once per matching
molecule/target/model configuration. Such a score cannot distinguish conformers
of that molecule and is not independent conformer-block evidence.

## Combining scores and allocating computation

Keep 3D, PLANTS and N-E columns separate. Orient/normalize against a shared,
target-specific control or frozen pilot reference, not separately inside each
block. Correlated models within a workflow do not each get a full independent vote.
With labels, calibrate on training controls and test held-out scaffolds; without
labels, rank consensus remains a prioritization heuristic, not binding probability.

For receptor aggregation, report best supported state AND state robustness.
Do not require success in every state: real ligands can favor one state. Also do
not reward a ligand merely for trying more states/poses; all candidates use equal
budgets and the same controls calibrate best-of-ensemble false positives.
Never assume equal crystal-state counts are Boltzmann populations or interpret
a raw minimum docking score as ensemble binding free energy.

Aggregate poses/conformers to unique molecules before estimating per-block tails.
Use a fixed upper-tail fraction (for example preregistered 5% means k=5 at 100,
k=10 at 200), report actual sample sizes, valid-pose rate, failures, disagreement,
state preference, geometry migration and time. Resample molecules for uncertainty;
preserve shared-molecule dependence in cross-block comparisons. Report adaptive
follow-up separately from the initial representative panel.

Expand high-priority, uncertain and disagreement blocks. Maintain a prespecified
random exploration allocation for low-ranked blocks. No good result among 100/200
samples is not proof a block lacks rare useful molecules. Time a small actual
panel first: cost scales with blocks, sampled molecules, conformers, microstates,
receptor states, seeds and retained poses.

## Validate the bias concern directly

1. Redock known ligands to their own receptors and cross-dock them into other
   states. Report top-ranked and best-returned pose recovery separately.
2. Exclude a test ligand's own complex, close chemical-family complexes and
   near-duplicate structures when constructing its evaluation receptor panel.
   This is harder than a diagonal-only redocking benchmark.
3. Compare single reference, diverse co-crystals and receptor-state ensemble
   strategies on held-out scaffolds and peptide/macrocycle controls where available.
4. Compare initial 3D-only filtering with unfiltered representative block docking.
   Record chemistry families lost before docking, not just survivors' scores.
5. Evaluate pose chemistry/stereo, ring closure, clashes, common receptor frame,
   known active/negative controls, repeated seeds and actual resource use.

A receptor ensemble reduces one source of template bias but does not guarantee
coverage of unseen induced-fit states. Add restrained receptor relaxation or
carefully validated MD-derived states only where cross-docking reveals a gap,
with an explicit separate protocol and cost budget.

## Current implementation boundary

Present: validated conformer blocks, nested sampler, 3D/multi-template workflow,
experimental structure diversity/alignment tools, PLANTS executable discovery and
path configuration. The current source explicitly marks PLANTS execution pending;
the existing execution adapter is Glide-specific. N-E is user-provided and its
callable interface has not been found in this repository. Pending: block-panel
export, PLANTS execution/results, receptor-state manifest, N-E adapter and joint
scorecards. Existing full-library blocking remains usable and unchanged.
