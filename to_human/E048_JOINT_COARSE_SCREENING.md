# E048 joint coarse screening

The former anchor-only filter could admit almost every conformer. The new optional
`coarse_constraints` policy requires all three explicit whole-ligand predicates,
then applies the selected anchor conditions before any pose generation. Every
library ID participates; there is no Top-K or forced rejection percentage.

| Predicate | Definition | Evaluated before seeds |
|---|---|---|
| Heavy atom ratio | Candidate count / crystal ligand count, inclusive interval | Yes |
| Shape extent distance | Euclidean difference of sorted principal RMS extents, divided by query extent norm | Yes |
| Chemical feature coverage | Sum over types of min(candidate count, query count), divided by total query feature count | Yes |
| Selected anchors | Type, distinct-assignment and compatible pair-distance/direction necessary conditions | Yes |
| Seed feasibility | Optimistic selected-feature score on existing rigid seeds | No: seeds required |
| Gaussian and E031 | Original winner selection and feature assignments on survivors | No |

The first three are NEW eligibility criteria. They are not rigorous upper bounds
on Gaussian similarity, detailed shape overlap, or protein binding. A molecule
that passed the former anchor-only rule may intentionally fail the new rule.
Eigenvalue extents omit detailed surface shape and chirality. Coverage counts
typed query features, not independent chemical interactions or activity.

No numerical defaults are silently imposed. A boundary guard of 1e-7 is applied
to normalized extent distance. Policies, masks and rejection reasons are saved.
Old anchor-only selections remain anchor-only: create a new selection preview.

## Pocket information

Use the existing classification to select explicit crystal-derived feature IDs:
CYS379 polar anchors; PHE433 aromatic/hydrophobic features; ASN376 polar features;
and other supported contacts with geometry evidence. ALL is within one query and
one stored candidate pose; do not combine independent candidate poses to satisfy
the rule. Shared-feature aliases do not become independent contact evidence.

Water bridges are optional hypotheses, not mandatory by default. Receptor atomic
clashes, pocket occupancy, protonation-dependent electrostatics and validated
candidate-protein contacts are NOT added by this change. A ligand-centered shape
filter cannot establish them. These require receptor preparation and placed-pose
checks; AF3 confidence does not substitute for pocket validation.

## Workstation procedure

Pull main and restart the chat backend with its existing profile and
`--allow-compute`. Reuse completed classified evidence in the same conversation.
Do not resume a code-hashed run created by an older implementation.

The following is an EXPLORATORY example, not an accepted cutoff set. It deliberately
combines three polar anchors and a crystal-derived aromatic feature. Adjust explicit
criteria based on inspected chemistry and measured retention, not a removal quota.

```text
Create a new selection preview from the completed QT9 classification in this
conversation. Require ALL of these exact anchors in the same pose:
8BJU:QT9:A:601/A0:HBD:30
8BJU:QT9:A:601/A1:HBA:29
8BJU:QT9:A:601/A2:HBA:39
8BJU:QT9:A:601/F12:pi_stacking:A:433:CG/CD1/CE1/CZ/CE2/CD2
Use minimum_score 0.5. Add coarse_constraints with heavy_atom_ratio [0.7, 1.3],
maximum_extent_distance 0.35, and minimum_feature_coverage 0.7.
These are exploratory eligibility criteria. Use no molecule cap.
Do not launch the full-library funnel yet.
```

After the preview completes:

```text
/coarse
```

This resolves the latest completed selection in the current conversation, samples
10,000 spread global IDs plus available saved positive/boundary rows, and does NO
seed generation or Gaussian evaluation. It reports retained/rejected conformers,
sequential rejection reasons, rejection fraction and whether the sampled 90%
target was met. Sampling is exploratory and deliberately supplemented; its rate
is not a statistically unbiased full-library estimate. The audit uses one CPU
process and excludes a full integrity scan; its latency is not full-run latency.

If selectivity is useful, run `/benchmark` for reference-versus-optimized pose
equivalence with the SAME rule, including current positive hits. Zero-positive
panels do not establish positive retention. Only then use `/funnel` for uncapped
whole-library counting and refinement. If coarse rejection remains weak, inspect
which requested features are discriminative; do not conceal failures with Top-K.

The low-level audit accepts a selection report path:

```bash
python -m aidd_agent.funnel_benchmark \
  --selection /absolute/path/to/selection/report.json \
  --output /absolute/path/to/new-coarse-audit \
  --count 10000 --coarse-only
```

## Evidence and limitations

User-reported E047: 288 rows; anchor invariant retained 288, seed feasibility
retained 199; NumPy22 median 1.526214 s versus 1.580669 s without seed feasibility.
No current reference positives; scalability gate false. This did not solve
pre-seed pruning. E048 workstation selectivity and full-library speedup remain
unmeasured. Neither 90% rejection nor RTX5090 acceleration is claimed.
