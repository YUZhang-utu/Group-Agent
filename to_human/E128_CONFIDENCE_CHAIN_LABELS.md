# E128: repeated AF3 confidence chain labels

The user supplied a completed live AF3 report for CIP2A and one ligand, with a
saved model CIF, explicit seeds 47..51, ipTM 0.84 and pTM 0.63. This is user-returned
execution evidence, not independent inspection of the Linux outputs or pose.

Repeated A/B identifiers do not imply thousands of chains or failed inference.
Full confidence arrays identify each token/atom, so repetition is normal there.
Official summary `chain_ids` instead has one entry per chain. The reported summary
contains repeated labels against two-element vectors and 2x2 matrices. The local
importer copies the source summary unchanged; the installed AF3 version and exact
source of this discrepancy remain unknown.

Changes:
- Confidence analysis validates unique labels and matching chain-level dimensions.
- `/results` compacts repeated labels into counts and warns about unverified matrix mapping.
- `task_report` adds read-time mapping diagnostics, including for existing results.
- New runs include confidence warnings and compact confidence-analysis artifacts.
- Raw reports, summaries, metrics, hashes and model outputs remain unchanged.

After pulling `feature/structure-guided-chat`, restart the Agent using the same
workspace, runtime and provider settings. In the existing conversation, use
`/results TASK_ID` with the completed task ID. No resume or AF3 rerun is needed
for these read-time diagnostics. Existing saved confidence-analysis files retain
their original contents; read-time structure confidence review uses the new checks.

95 targeted software tests passed. No new live-provider or GPU calculation was
performed. The mapping warning concerns chain labels, not proof that recorded
scalar scores or the predicted coordinates are wrong. Do not infer affinity from
these scores or silently deduplicate an unverified array into matrix labels.

References: [AF3 output format](https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md),
[upstream summary writer](https://github.com/google-deepmind/alphafold3/blob/main/src/alphafold3/model/confidence_types.py).
