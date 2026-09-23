# E055: score contracts and workstation validation

## What changed

New recommendations use `optional_normalization: fixed_budget`. The coordinator
computes `optional_budget` once from the initial sum of independent weights and
family weights (one weight per family, not per member). Edits retain this budget.
For the supplied 20-mode proposal the budget is 9.75.

```
contact = (sum(w_i * match_i) + sum(w_g * max(group_matches))) / optional_budget
composite = gaussian_weight * gaussian + optional_weight * contact
            + occupancy_weight * occupancy - soft_exclusion_penalty
```

Fixed-budget dimension coefficients must sum to one. Weight totals exceeding the
budget are rejected, never clipped or automatically renormalized. A zero-weight
initial proposal receives budget 1. Legacy designs without these fields preserve
their weighted mean; there is no silent migration of old sealed designs.

A budget is a scale choice, not an independently validated theoretical maximum.
Separate new recommendations may start with different budgets; use the SAME
explicit budget when comparing variants. Deleting a contact leaves unused budget,
so the achievable score can decrease. Deliberately increasing the budget changes
the scale and must be reviewed. Deleting anchors may also change the existing
at-least-one-anchor eligibility rule. Re-adoption with a null minimum_pose_score
recomputes the exploratory crystal-derived threshold; score invariance alone is
not proof of invariant library membership. Freeze a documented threshold for a
controlled ranking comparison, and separately report any recalibrated threshold.

Each spatial group optionally accepts `ambiguity_margin` in Angstrom, falling
back to global `spatial_ambiguity`. Radius controls eligibility. Among eligible
regions, the nearest wins only if its distance gap to EVERY other eligible region
exceeds the larger of their two margins. Otherwise the atom is ambiguous and
counts toward neither region. This conservative pair rule is order-independent.
No automatic variance-to-radius estimator or target-specific constants are added.

New consensus anchors, adviser context (including older surveys), and adopted
anchor details expose `protein_part`, `ligand_role`, and `interaction_type`.
Canonical protein atom names distinguish backbone O from sidechain CG1/CG2;
HBD/HBA refer to ligand donor/acceptor. Evidence IDs and coordinates are unchanged.

## Multi-template audit

The existing consensus funnel scans each template separately, unions molecule
membership, and retains actual poses per molecule/template/anchor mask/spatial
signature. It never averages template coordinates. There is no calibrated global
cross-template ranking. Within each template, representative poses compete on
their composite score. Crystal threshold derivation uses the largest same-pose
Gaussian score across selected templates, without statistical calibration.
Reports now expose this policy explicitly. Do not interpret raw template scores
as comparable affinities or treat a union as a calibrated best-template ranking.

The local supplied-proposal replay evaluated 45 crystal identity poses against 11
templates (495 scores). Nonself medians ranged from 0.2228 (5HMK, 50 heavy atoms)
to 0.3569 (3JZK, 31 heavy atoms). Thus template distributions differ, but this
diagnostic does not support a universal larger-template advantage. Self-template
pairs are excluded from distributions and winners; exact winner ties split credit.
Related compounds remain in the panel: this is not held-out or library calibration.

## LYS94 evidence review

Mode `pocket:a6190c97bee1efc7` has five structures and three chemotypes. Recorded
ionizable-feature-center to protein NZ distances are:

| Complex | Feature-center distance (Angstrom) | Nearest ledger oxygen distance (Angstrom) |
| --- | ---: | ---: |
| 4ERE:0R2:A:201 | 2.965 | 2.662 |
| 4ERF:0R3:A:201 | 3.868 | 3.182 |
| 4OAS:2SW:A:201 | 4.260 | 3.669 |
| 4QO4:35S:A:201 | 3.332 | 2.733 |
| 4QOC:35T:A:201 | 3.425 | 2.716 |

These are different distance definitions, not contradictory measurements. The
nearby ligand oxygen atoms carry RDKit NegIonizable features; their CCD-perceived
formal charges in this ledger are zero. Salt-bridge labels therefore express an
ionization hypothesis, not validated protonation or energetic necessity. Solvent
exposure, binding contribution, and Asp/Glu library enrichment remain unmeasured.
The evidence does not justify treating 0.80 as a measured energetic weight.

## Run on the workstation

Stop Chat, update the existing branch, and preserve the same storage root:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
conda activate aidd-workstation
git switch feature/structure-guided-chat
git pull --ff-only origin feature/structure-guided-chat
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python -m pytest tests/test_score_contract.py tests/test_contact_groups.py tests/test_consensus.py -q
python scripts/check_english.py
```

Restart Chat with its original command, storage, runtime and provider settings.
Completed consensus can be reused; issue a NEW `/recommend` rather than resuming
an old code-hashed task. Verify `optional_normalization: fixed_budget` and a
positive `optional_budget` in its returned recommendation. Do not expect exactly
9.75 if the model selects different contacts or weights.

For a controlled comparison of the supplied proposal, keep its contacts/templates
and explicitly edit `optional_normalization` to `fixed_budget` and `optional_budget`
to 9.75. Subsequent weight edits must retain 9.75. Chat design edits directly create
an adopted design with crystal checks; `/adopt` afterward would accept the original
proposal instead of preserving the edit. Adoption is not a library scan.

For read-only diagnostics, use the actual recommendation task's saved report.json,
not a pasted UI JSON (the latter omits sealed sources and the survey path):

```bash
python -m aidd_agent.score_audit \
  --recommendation /absolute/path/to/recommendation/report.json \
  --output /absolute/path/to/new/e055-audit
```

Replace both paths. Output must not already exist. This reads original evidence
and writes diagnostics only; it does not adopt, alter a source, call an LLM, or scan
the library. Check `fixed_budget_water_removal.passed`, `template_distributions`,
`selected_contact_evidence` (including atom pairs and chemistry assumptions), and
`pose_scores`. A passed score audit does not mean a funnel or activity test passed.

## Local verification and remaining scientific work

- Full regression: 426 passed, 2 skipped before the final audit-output enrichment;
  focused final checks are recorded in research_log.md.
- Fixed-budget water-removal invariance: passed for all 45 local crystal poses.
- Explicit fixed-budget adoption fixture: 45 crystal controls, no selected-template
  failures; production user design was not adopted.
- Fixture source: data/e055-user-proposal-fixture/report.json.
- Full local diagnostic: data/e055-score-audit-sealed/report.json.
- Adoption fixture: data/e055-fixed-adoption-fixture/report.json.

Next scientific validation needs reviewed spatial atom selections, NH/library
retention, independent positive controls and a representative background panel.
Compare dimension ablations on the same panel and score budget before choosing
0.4/0.4/0.2 or any other mix. Pocket occupancy is a plausible signal, not yet shown
to be the strongest. Cross-template calibration and per-region uncertainty
estimation remain explicitly unimplemented, rather than fitted to these crystals.
