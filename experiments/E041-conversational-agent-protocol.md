# E041: persistent conversational workstation agent

Goal: a local English chat UI with selectable model provider, persistent sessions,
background scientific tasks, status/results queries, cancellation and explicit resume.
Reuse the owned Project planner/executor and scientific receipts. User requests may
be in any language; maintained UI and model replies remain English.

Acceptance: conversation clarification retains context; query-only messages never
launch compute; task IDs cannot cross sessions; jobs survive UI reloads; crashes
mark in-flight tasks interrupted without silent rerun; status/results come from
actual reports; task queue serializes heavy jobs; cancellation stops subprocess
groups; resume uses the original sealed plan; endpoints require a local bearer
secret and reject cross-origin mutations. Test with fixtures and a real local HTTP
server; real models/GPU remain workstation acceptance. Record stage capability
boundaries, including candidate retrieval versus docking and biological validation.

Do not silently map new targets to WEE1. Generic protein/AF3 workflows are already
supported, while calibrated full-library queries remain WEE1-only. Existing CLI
pre-docking QC is not an independently validated docking engine. Missing workflow
inputs/capabilities must be surfaced, not simulated as completed scientific work.
