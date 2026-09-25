# Automatic structure-to-library continuation

E074 adds a persistent coordinator around the existing asynchronous Chat jobs.
It starts from an existing task, which may still be queued, and advances completed
diversity -> consensus -> recommendation -> budget stages. Choose a final goal of
recommendation or budget. It does not run new scientific algorithms or docking.

## Update and start

Update `feature/structure-guided-chat`, restart Chat with the same storage/runtime/
provider settings, and refresh the page. Use the session holding the source task.
For budget execution the server must retain `--allow-compute`. Do not change an
active sealed search's code/runtime and assume that job can resume unchanged.

Preparation-only automatic continuation:

> Automatically continue from diversity task TASK_ID through the aligned consensus and consensus recommendation. Use reference ligand 6Q9L:HTZ:A:201 and target chain A. Stop after the recommendation; do not screen the library.

Only use these reference identifiers if they are the intended, reported MDM2
instance. For another target supply its actual identifiers.

Explicit screening and delivery authorization:

> Automatically continue from diversity task TASK_ID through consensus, recommendation and budget screening. Use reference ligand 6Q9L:HTZ:A:201 and target chain A. Retrieve 1000000 unique molecules and export the first 100000 with their original names and original and posed MOL2 files. Keep the configured workers and chunk size. Stop after this first delivery.

Starting from an existing recommendation or design skips prior stages:

> Automatically continue from recommendation task TASK_ID to budget screening and the first 100000-molecule delivery, retrieving 1000000 unique molecules. Reuse the existing recommendation and stop after delivery.

For a new target, first request its protein/structure-diversity task as in E073.
The agent can register a requested chain against that newly queued task ID. A chain
does not need to wait in the browser; the Chat server must remain running.

## Monitor and control

> Show my automatic structure workflow chains, current task IDs and any blocked dependencies. Do not start new work.

> Pause automatic continuation for chain CHAIN_ID. Leave its current task running.

> Resume automatic continuation for chain CHAIN_ID.

> Cancel future continuation for chain CHAIN_ID. Leave existing results intact.

The task panel has a Structure workflow card with the chain ID, goal, current task,
source/child history and errors. State persists in `chat/conversations.sqlite3`.
Paused chains remain paused after restart. Completed dependencies can advance after
restart; interrupted scientific tasks are not silently resumed. Resolve and resume
the scientific task separately, then resume a blocked chain.

## Failure semantics

- Missing reference identity or failed preparation blocks continuation. No default
  ligand, chain, score cutoff or target is invented.
- Pausing/cancelling the chain affects future stages only. Cancel the current task
  separately if that is also wanted.
- A crash between dispatch and saving its receipt marks the chain `uncertain`.
  It cannot automatically resume. Inspect task receipts first; cancel the old chain
  and create a new one from the confirmed task after reconciliation. Do not restart
  from the original source blindly.
- A blocked chain does not automatically retry a failed scientific task. Existing
  seals, source integrity and compute policies still govern that task.
- Budget completion means the adapter completed, not that the requested molecule
  count was necessarily available. Inspect shortfall, review flags and validation.
- Budget retrieval queries the library index and scores retrieved candidates;
  this is not exhaustive pose scoring or validated biological enrichment.

## Workstation acceptance

First test the recommendation-only goal using an existing diversity task. Verify
one consensus job, then one recommendation job, and no search job. Pause before a
dependency finishes and check that the next stage is withheld. Restart Chat and
confirm the chain ID/state survive. For a budget goal reuse an appropriate source
and explicitly authorize compute. Do not rerun a completed full-library search
solely to test this scheduler. Local fixture tests do not establish live provider,
PDB network, desktop or production-search success.
