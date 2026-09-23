# E056: contact evidence, self-matching and assignment diagnostics

## Scope and changes

New consensus observations retain source instance, native feature index, ligand
atom indices, extraction method, criterion version and unvalidated assumptions.
Representative atoms are always namespaced by their representative query/source
anchor. Existing consensus IDs do not change merely because metadata was added.
The complete atom ledger remains separate from sparse feature extraction.

Recommendation and adoption reports now expose a scoring_contract, including the
actual matching path, denominator, dimension coefficients, and post-assignment
family-max semantics. The Chat summary displays this contract. Old direct scoring
interfaces retain their legacy semantics; there is no claim of cross-target score
calibration. Scientific weights, sigma and angular power are unchanged.

## Three diagnostic layers

1. Feature identity: every ligand feature matches the same query under the identity
   transform. Valid features should score approximately one.
2. Source observation: each recorded source contact matches its original ligand
   feature alone in the aligned frame. This isolates extraction, transformation,
   direction and mapping consistency from multi-anchor assignment competition.
3. Selected consensus: each crystal ligand matches the selected pooled consensus
   modes with the actual one-to-one feature assignment. Missing modes and imperfect
   matches are expected and are not classified as implementation failures.

Reports contain feature type names, direction kinds, unmatched counts, score
quantiles, assigned spatial/angular factors, raw rows and per-type weighted
contributions. Unmatched scores are zero, while unavailable geometry factors are
null. Geometry-factor summaries state their denominator and exclude null values;
score summaries include unmatched zeros. Family contribution ties split credit.
Signed and axial angles are feature-direction comparisons, not direct validation
of candidate hydrogen bonds or ring interactions with protein atoms.

Self-match success does not prove that angular penalties on new candidates predict
activity. The controlled opposite-direction regression distinguishes signed from
axial behavior without choosing new empirical compensation weights.

## Suspected shared evidence

The audit compares only observations from the SAME source query instance. A pair
is flagged when atom-index sets overlap and aligned feature centers lie within
the supplied distance. The default 0.5 Angstrom is an exploratory review threshold,
not an automatic scientific merge rule. Report both protein-partner and interaction
type agreement. Missing verified indices are unresolved, never negative evidence.

The output is a pairwise graph, not connected-component scoring groups. A-B and
B-C relations do not merge A with C. Different queries never share atom identities
just because their local atom numbers happen to match. Different protein partners
and interaction types can share a ligand group without being duplicate physics.

Old surveys can be audited without rebuilding: mappings are recovered only from
source-sealed native manifests and their checksum-verified query packages. Missing
packages are reported as unresolved. Recorded new mappings must agree with the
sealed source. Sources are never modified by the audit.

## Optional-group assignment comparison

Production still performs individual-feature assignment followed by family max.
The diagnostic collapses each optional family into a row whose score for a candidate
feature is the best member score, then solves the weighted one-to-one assignment
between units and candidate features. Independent contacts are singleton units.
Small exhaustive tests establish the optimum of THIS optional-only objective.

This alternative admits one contributing candidate per optional unit; it is not
an automatic replacement for the production hard-rule, hit-mask or multiple-contact
semantics. Hard-rule designs are explicitly marked not applicable. The report
shows both objective values and the alternative member/feature selections. Feature
uniqueness is not atom uniqueness because one atom can generate multiple features.

## Workstation commands

Stop Chat, pull the existing branch and activate the environment:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
conda activate aidd-workstation
git switch feature/structure-guided-chat
git pull --ff-only origin feature/structure-guided-chat
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python -m pytest tests/test_contact_diagnostics.py tests/test_contact_groups.py tests/test_score_contract.py tests/test_consensus.py tests/test_prompt_workflow.py -q
python scripts/check_english.py
```

Run the evidence audit against the saved recommendation task report, not copied
chat JSON. Replace both absolute paths; output must be new:

```bash
python -m aidd_agent.evidence_audit \
  --recommendation /absolute/path/to/recommendation/report.json \
  --output /absolute/path/to/new/e056-audit \
  --duplicate-distance 0.5
```

Open `review.md` for readable distributions and `report.json` for source-level
evidence, score factors, contributions, pair relations and assignment diagnostics.
Require `feature_identity_passed: true`, `source_recovery_passed: true` and zero
unresolved source mappings for full source-recovery acceptance. Null is not pass.
This command performs no LLM request, adoption or library scan.

The existing E055 score audit remains available for template distributions and
fixed-budget deletion invariance. Reuse the same saved proposal to compare runs;
a new LLM recommendation may change anchors and weights and is not a controlled
before/after comparison. Do not resume old code-hashed execution tasks.

Restart Chat with the SAME original workspace and runtime configuration. Existing
consensus supports the audit recovery path. Rebuild consensus only when permanent
new observation metadata is desired, then create a fresh recommendation tied to
that report. Do not mix proposal and survey sources manually in production.

## Scientific boundaries

Do not automatically compensate directionality, merge suspected pairs, alter
production assignment, or interpret crystal coverage as active/decoy enrichment.
Review the emitted conflicting cases first. Biological recall testing remains
deferred at the user's request. Profiling is local crystal matching, not end-to-end
library throughput or a benchmark of larger target-specific anchor sets.

## Local results

Full suite: 439 passed, 2 skipped. Final focused command above: 92 passed.
An offline raw-structure rebuild retained 45 complexes, 119 modes and unchanged
anchor IDs, with zero preparation failures. The real audit recovered all mappings:
878 identity features and 283 source-observation rows scored approximately one.
The selected-consensus layer had 900 rows with 456 unmatched; these are not expected
to all match. There were 63 suspected shared-evidence pairs, and the alternative
optional objective improved 13/45 crystal poses (maximum raw numerator gain
0.0828388993). Production assignment is unchanged. See E056_VALIDATION_SUMMARY.json.
