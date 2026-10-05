# E098: MEDCHEM Chat and completed CLI docking adoption

Engineering objective: retain completed E097 results and control their inspection
and analysis through the existing Chat, without launching PLANTS again.

Protocol before implementation: bind an attachment to the active Project and session,
record the immutable report digest, inspect bounded counters without loading poses,
and queue explicit analysis through the existing task executor. Analysis must verify
all original output seals and reject changed reports or paths outside the Project.
Repeated attachment/analysis must reuse existing IDs. Test partial panels, changed
inputs, cross-project paths, provider tool dispatch and no engine invocation.

UI: use the supplied med.png unchanged, a teal molecular masthead, locally rendered
illustrative cyclic-peptide 2D structures, and hover/focus/click workflow coverage.
Keep existing session, provider, viewer and task controls. Check desktop/mobile
layout, keyboard access, HTTP assets and existing regression tests.

The user reports full docking completion on 2026-10-05. Final workstation path and
receipt have not yet been inspected. No new scientific experiment, affinity result,
block rejection policy or live-provider acceptance is implied by local tests.

Outcome (2026-10-05): the user supplied the final CLI receipt: 8,445 attempted
jobs of 8,445, zero failures, 821,619 scores, first-job gate passed. This is
user-reported workstation evidence; artifacts were not transferred here.
Local regression: 743 passed, 6 optional-dependency skips. Nine adoption tests
cover a real background Python analysis worker and a synthetic provider tool loop;
PLANTS calls occur only in fixture construction through a mock. Desktop/mobile
browser checks passed and screenshots were visually inspected. These are
confirmatory engineering checks, not a scientific comparison of partitions.
