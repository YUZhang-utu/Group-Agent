# E074: durable structure workflow continuation

Hypothesis: a persisted coordinator can advance already-authorized structure stages
without repeated user task-ID handoffs, while refusing duplicate uncertain dispatch.

Protocol: fixture a diversity task; run coordinator ticks through consensus,
recommendation and budget job receipts. Repeated ticks must wait for queued tasks.
Test pause/restart, failed dependencies, compute disabled, session ownership,
duplicate registration and loss of dispatch receipt. Scientific execution is mocked;
task storage, context reads and scheduler transitions use the production code.
Run the full regression because the shared worker loop is modified.

Acceptance: exactly one dispatch per completed stage; no progress past failed or
ambiguous input; uncertain dispatch never retries automatically; completed goals
stop. Existing jobs remain separately resumable/cancellable. No real library run.

Focused confirmatory result: 25 tests passed. Live provider interpretation and
workstation service lifecycle remain acceptance checks, not established by fixtures.
