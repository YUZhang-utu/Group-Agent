# E090 matched-target findings

Input D:/agent/pair_metrics.csv matches the supplied SHA256 exactly:
5167b31a3896c342b6ec1dceadc2e668474e5812e2327f636b747ed16cfb297d.
All 7,176 pair rows were read; 2,107 unique conformers were recovered. No
same-molecule pairs occur. This validates the supplied table identity, not a
fresh recomputation of original coordinates (profiles/MOL2 were not supplied).

## Matched comparison

For each query/target/rotation, take the median over reference molecules and
subtract the median reference-reference difference within that same target.
Then take the median over candidate rotations/targets per query, and summarize
100 equally weighted queries per band. All observed rotations are retained.

| Band | Property delta | Steric delta | Local steric maximum delta | Backbone RMSD delta (A) |
|---|---:|---:|---:|---:|
| <=35 | +0.00205 | -0.00056 | -0.00523 | +0.02048 |
| 35-45 | +0.00313 | +0.00451 | -0.00381 | +0.02643 |
| >45 | +0.00395 | +0.00311 | -0.03391 | +0.00603 |

These medians show no large central shift relative to the selected class
controls. They do not prove equivalence, identify a universal threshold, or
show that every candidate is acceptable. The p90 query-level backbone excess
is +0.256 A, +0.308 A and +0.475 A respectively: tails still need review.
Property and steric metrics are scaled descriptor quantities, not angstroms.

The original pooled <=35 RMSD of 1.204 A versus control 0.950 A mixes different
targets and weights. Comparing the 123 target classes shared between <=35 and
35-45, the median target-paired RMSD difference (35-45 minus <=35) is -0.00033 A.
The pooled values do not establish that <=35 is worse.

There are only 42 targets represented in all three bands. Within this common
subset, median target-specific backbone excesses are -0.069 A (<=35), -0.050 A
(35-45), and +0.232 A (>45). This exploratory subset suggests that the nearly
zero overall >45 delta should not be interpreted as evidence for unrestricted
relaxation. It is a subset comparison, not a significance or population claim.

## Candidate-specific fit, not band-wide acceptance

For a deliberately descriptive check, compare each candidate/reference median
to its own target's observed control maximum on all four metrics simultaneously.

| Band | Queries with at least one candidate inside this envelope | All matched candidates inside |
|---|---:|---:|
| <=35 | 93/100 | 0/100 |
| 35-45 | 89/100 | 4/100 |
| >45 | 87/100 | 1/100 |

This supports investigating candidate-specific placement, not assigning all
angle-compatible target classes. These are NOT acceptance or recall rates.
Multiple target/rotation opportunities increase the chance of finding a match.
The envelope itself is weakly estimated: 599 targets have just three control
pairs and two have only one. Six targets have no controls; their six candidate
rotations are excluded from matched analysis and exported separately. The
"all" column refers only to matched candidates, not those without controls.

## Implications for the next block revision

1. Keep definite cis/trans and handedness as compatibility constraints.
2. Prioritize <=35 candidates and test 35-45 as an extension. Do not use the
   current aggregate results to automatically move all >45 records.
3. Select the specific target and directed rotation using backbone, continuous
   chemical properties and local steric checks jointly. Exact sequence identity
   is unnecessary, but hydrophobic annotation alone is insufficient.
4. Calibrate against additional diverse references and held-out definite
   members before using any observed class envelope as a routing cutoff.
   These first-indexed three-molecule reference sets are too small to define
   every class's limits. Some backbone classes may themselves be broad.
5. Preserve uncertain/rejected candidates in the exhaustive special pool.
   Measure its actual search time and held-out retrieval before claiming final
   >95% validated block coverage. Do not multiply this class-balanced panel's
   fractions by the full-library conformer population.

No boundary membership has been changed by this analysis.

## Reproduce and inspect

```bash
python scripts/analyze_e089_pairs.py --source pair_metrics.csv \
  --output NEW_ANALYSIS \
  --expected-hash 5167b31a3896c342b6ec1dceadc2e668474e5812e2327f636b747ed16cfb297d
```

The analysis exports matched_candidates.csv, query_balanced.csv,
missing_controls.csv, common_targets.csv and band_target_comparisons.csv,
plus the machine-readable report. The delivery ZIP also contains the original
pair table and a SHA256 manifest. All analyses are exploratory descriptions of
the supplied workstation panel; they are not binding affinity or docking tests.
