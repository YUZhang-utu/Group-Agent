# E128: AF3 confidence chain labels

## Protocol

The user supplied a live completed AF3 report with many repeated A/B entries in
summary `chain_ids`, but two-element chain scores and 2x2 pair matrices. The local
adapter copies the source summary unchanged. We have not inspected the user's
Linux files or the installed AF3 version.

Hypothesis: the report needs chain-label validation and compact presentation;
repeated labels alone do not show inference failure. Official AF3 output documentation
defines summary `chain_ids` per chain and full `token_chain_ids`/`atom_chain_ids`
per token/atom. Do not guess which repeated array was supplied or silently deduplicate
it into matrix labels.

Preserve raw files, numeric metrics, and hashes. Add a read-time assessment of
chain-label uniqueness and chain-level dimensions, compact repeated arrays with
explicit counts, and warn that matrix mapping is unverified. Apply read-time checks
to existing task results as well as new confidence explanations. No GPU rerun.

Validate synthetic repeated/unique/missing/mismatched labels, preserve input data,
and run the existing prediction, structure-chat, domain-tools, and reliability
tests plus the English content check. These are software checks, not independent
validation of the user's structure.

Reference: https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md

## Results

Confirmatory software checks: 95 tests passed across confidence mapping,
prediction import, structure Chat, domain tools, task reliability and prompt
workflow suites. The initial integration fixture lacked its required plan path;
correcting the fixture allowed the existing read-only structure review to run.
Raw summary bytes, hashes, metrics and completed status remain unchanged.

Current upstream `model/confidence_types.py` derives summary chain_ids by unique
first occurrence of token labels. The user's installed code/version and raw
files have not been inspected, so the origin and granularity of the repeated
summary array remain unverified. Count summaries do not assign matrix order.
