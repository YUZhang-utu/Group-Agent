# E104: recorded block scores to conformer search and docking

## Protocol (2026-10-06, before execution)

Engineering integration, followed by exploratory workstation evaluation. Reuse
the completed 821,619-pair EquiScore panel; do not repeat inference or docking.
Adopt the sealed analysis and its ChemPLP baseline into the owning conversation.
Select Top-5 or Top-10 blocks within one scheme/receptor, using EquiScore,
ChemPLP, union or intersection. The block statistic remains the mean of the
best 10 distinct molecules (best conformer per molecule).

Use the existing confirmed MDM2 consensus query. Search only conformers owned
by selected blocks. The delivery budget is 100,000 conformers, explicitly
allowing multiple conformers of one molecule. Count unique source conformers,
unique molecules and conformer/receptor pairs separately. Never expand a
selected molecule into conformers outside the selected blocks. Report shortfall
when the qualified pool is smaller; never fill it with unrelated candidates.

Reuse contact-first 3D refinement and the original source chemistry. Preserve
the original query, source-to-artifact identity mapping, selected memberships,
scores, transforms and downstream PLANTS pose references. Dispatch through
the existing sealed Project task runner, with idempotent requests and receipts.
The explicit search-and-dock request authorizes both steps.

## Verification

Fixtures must cover ownership, changed seals, conformer versus molecule counts,
model union/intersection, deterministic ties, membership restriction, source
mapping failures, missing/invalid search poses and repeated task submission.
Local tests do not establish scientific recall, affinity or workstation success.
Real acceptance requires importing the workstation reports, choosing the owned
MDM2 query task, and completing a new search/docking run in the Agent.

## Local implementation outcome

Implemented owned, asynchronous adoption of existing full EquiScore analysis,
ChemPLP baseline inspection, union/intersection selection, source-hash mapping to
search artifact IDs, conformer-preserving contact-first ranking/export, original
reviewed PLANTS configuration inheritance and saved per-conformer results.
Chat tools, explicit JSON fallback commands and drawer controls share contracts.
Repeated requests reuse receipts, including failed submissions; imported tables
do not trigger another model run. Missing artifact mappings fail before search.

Fourteen new engineering tests passed, including a synthetic six-conformer /
three-molecule search-to-PLANTS continuation, explicit 99,994 shortfall against a
100,000 request, repeated receipt reuse, changed source rejection and prompt
submission stopping after its queued receipt. This is confirmatory evidence for
the tested software contracts only. No real 100,000-conformer run was launched.

The real subprocess compute-gate test exposed duplicate `Blocked` exception
identities under `python -m`; the exception now lives in a shared module so
missing prerequisites are reported as blocked, not failed. Import completion
is published only after both EquiScore and ChemPLP views are copied.

Final combined regression: 139 tests passed across campaign, block results,
block ranking/evaluation, EquiScore, Chat, domain-agent, prompt workflow and
budget-search suites. JavaScript syntax and the English-content guard passed.

## Configuration-block recovery protocol (2026-10-06)

User task f0656f4d2c584dc4 stopped before scientific execution: the runtime has
search.batch and sampling_profile but only receptor_tools_profile, not a separate
plants_profile. Test reuse of that trusted engine configuration with the sealed
original receptor/site settings. Engine hashes must still match. Add explicit,
idempotent retry_config for configuration-blocked campaigns only, preserving the
old receipt and creating a new sealed plan. Reject retries after scientific
artifacts exist, for other task types/statuses, or across conversations. Validate
with synthetic pipeline and queue tests before delivery; no remote run implied.

Configuration recovery validation: 66 related tests passed, including both
PLANTS configuration routes through synthetic search/export/docking, engine
mismatch rejection and idempotent retry with session/artifact/status guards.
English guard passed. Real workstation retry remains pending.

Follow-up protocol: the real executor creates campaign/blocks before checking
runtime, so the empty-directory retry guard was too strict. Permit only that
known empty directory; reject nested files, other entries and completed stage
receipts. Extend the retry regression to match the executor directory layout.
