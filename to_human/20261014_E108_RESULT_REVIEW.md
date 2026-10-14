# E108 returned-output review, 2026-10-14

## Evidence and integrity

User supplied six files under D:/agent/aidd: arm_summary.csv, audit.json,
pairwise_overlap.csv, recovery_work.svg, reference_stage_coverage.csv and
top5_top10_gains.csv. historical_recovery.csv and report.md were not supplied.
Do not rerun the scientific searches; request the missing analysis outputs.

The audit reports complete historical-reference analysis for twelve arms and one
current-query hash. Its implementation hash matches the local E108 source.
All 60 stage rows satisfy partition accounting and nested count inequalities;
million-reference export counts agree with arm_summary.csv. There are 132 pairwise
rows. These are local consistency checks on returned outputs, not an independent
revalidation of the remote raw candidate files or registry.

## Primary six full-budget arms

All export 100,000 conformers. Counts below are exact conformer matches, not
activity hits or molecule-level reference recovery.

| Scheme/scorer | Examined conformers | Historical Top-100,000 recovered | Historical Top-1,000,000 recovered |
|---|---:|---:|---:|
| E094/ChemPLP | 257302 | 366 (0.366%) | 3809 (0.3809%) |
| E094/EquiScore | 139603 | 294 (0.294%) | 2980 (0.2980%) |
| E095/ChemPLP | 203971 | 210 (0.210%) | 1699 (0.1699%) |
| E095/EquiScore | 197604 | 396 (0.396%) | 4085 (0.4085%) |
| E096/ChemPLP | 214306 | 147 (0.147%) | 1186 (0.1186%) |
| E096/EquiScore | 176353 | 439 (0.439%) | 5199 (0.5199%) |

E096/EquiScore leads the Top-100,000 and Top-1,000,000 historical recovery endpoints,
and has the highest million-reference overlap per examined conformer among these
six arms (about 2948 recovered per 100000 examined). This ratio describes the
observed operating point, not an interpolated equal-work experiment or runtime.
The ranking is not stable across reference cutoffs: E094/ChemPLP leads Top-10000
recovery (49 versus E096/EquiScore's 33); E095/ChemPLP leads Top-1000 (6 versus 3).
Low counts prohibit treating those differences as robust general superiority.

## Where historical references disappear

For E096/EquiScore/Top-10 and the million reference:

- 991666 lie outside the ten selected blocks.
- 8334 lie inside selected blocks; 2979 are in blocks not reached before stopping.
- 5355 lie inside the four actually visited blocks.
- 52 of those do not qualify under the current search; 5303 qualify.
- 104 qualified references are removed by output truncation; 5199 are exported.

Conditional on being in visited blocks, 5199/5355 = 97.09% are retained. This is
not overall recall: overall historical-reference recall remains 0.5199%.
Membership/visited-block coverage dominates this endpoint's loss, not final
truncation. However, the historical/reference universe and current regular-block
eligibility are not yet reconciled. References outside selected blocks must not
all be called recoverable routing mistakes; inspect their full membership first.

## Limits and diversity of selections

Four Top5/Top10 pairs have exactly identical conformer and molecule outputs:
E094/ChemPLP, E095/EquiScore, E096/ChemPLP and E096/EquiScore. They already stop
within five blocks. Only E094/EquiScore and E095/ChemPLP expand the search at Top10.
The former adds 22992 and loses 2081 conformers after reranking; the latter adds
29331 and loses 6264. A larger allowed block count does not imply set containment.

Among the fifteen pairs of six Top10 arms, fourteen have zero exact conformer
overlap. E095/ChemPLP versus E096/ChemPLP shares 2250 conformers. Molecule overlap
is also mostly zero or very small. This indicates strongly different selected
regions, not proof of complementary high-quality chemistry or an identity bug.
Inspect full block priorities and parent memberships before causal interpretation.

## Conclusions and next steps

Current results do not support a claim of high historical-candidate retention.
They also do not show that newly found candidates are inferior: historical ANN
quota-tier RRF and current threshold/RRF differ, with no activity truth available.
No random-control replay or matched actual-runtime evidence was supplied. All
Top10 times reuse cache; no cold speedup can be inferred.

Request historical_recovery.csv to inspect molecule-level reference recovery and
common-prefix sensitivity. Next reconcile the historical references against the
complete current eligible membership and compare the frozen query/filter policies.
Only then run a controlled complete-universe or fixed-sublibrary random-order
benchmark. Do not tune on this reference and describe the same data as held out.

## Follow-up: supplied recovery table and execution-policy audit

The user subsequently supplied historical_recovery.csv. All full-output conformer
hit counts match reference_stage_coverage.csv. Top10 molecule-level reference
Top-100000 recovery is E094/ChemPLP 694 (0.694%), E094/EquiScore 630 (0.630%),
E095/ChemPLP 389 (0.389%), E095/EquiScore 674 (0.674%), E096/ChemPLP 308 (0.308%),
and E096/EquiScore 944 (0.944%). Molecules use first occurrence in the historical
conformer ranking. Low overlap is therefore not explained solely by choosing
different conformers of the same molecules.

Code-path inspection identifies a larger difference than quota-tier ranking:

- budget_screen.retrieve performs whole-index USRCAT FAISS retrieval, per-template
  unique-molecule selection and RRF, then collects stored conformers for the
  selected molecules. This is the historical budget-search path recorded for
  PROMPT-59addb9ffd6444c3, not an exhaustive threshold scan.
- budget_screen.make_queries defaults to budget mode; threshold mode is explicitly
  set for current E106 arms. Both share preselection_full.compute but branch there.
- In budget mode, preselection_full.compute bypasses coarse/anchor-bound and
  mandatory-geometry rejection. pose_representatives does not enforce rule_passes
  or minimum_pose_score as rejection; the score threshold is an annotation.
  Pocket exclusions remain applied in both modes. Budget mode is not unscored or
  free of physical filtering.
- Current threshold mode enforces adopted coarse, geometry, anchor-rule and pose
  score conditions. Tests cover equivalence to the full-library threshold funnel,
  not equivalence to the historical ANN budget pipeline.
- Historical conformer export additionally uses the saved quota tier; current
  ranking uses the current searched pool without that tier. RRF depends on its
  population, so identical raw scores do not imply identical cross-pool ranks.

These are implementation and recorded-lineage findings. Actual old versus new
parameter values/code signatures still require the saved remote protocol/design
artifacts; current source alone does not authenticate historical execution.
Most observed reference loss occurs before those scoring differences can act:
only 5355 of the million-reference conformers are in E096/EquiScore's visited blocks.
Calling changed thresholds the main numerical cause would contradict stage counts.

The intended scientific comparison should hold the evaluator fixed and change only
routing. Freeze either the current threshold evaluator on a common complete test
universe, or replay the historical budget evaluator consistently; do not silently
change production policies or overwrite existing receipts. First inspect complete
reference-to-block membership and raw protocol differences. No new scientific run
or code change was launched by this review.
