# PDB to consensus to molecule-library screening

Update: E074 adds opt-in persistent automatic continuation. See
[the automatic chain guide](E074_AUTOMATIC_STRUCTURE_CHAIN.md). The stepwise
interface below remains available; the E073 scheduler limitation describes that release.

The domain agent now has a read-only `structure_workflow` tool that identifies
current-session stages, actual source-run links, proposed reference IDs, reference
ambiguities, similarity artifacts and supported next actions. Existing executors
still perform the scientific work. No scores, thresholds or retrieval methods were
changed by this integration.

## Scientific chain

| Stage | Existing implementation | What Chat should inspect |
|---|---|---|
| Target and experimental structures | Protein lookup and structure_diversity, with supplied reference PDB | Target identity, discovery coverage, quality/site exclusions |
| Ligand differences | Morgan radius-2, chirality-aware Tanimoto matrix; quality-seeded farthest-first proposal | ligand-similarity.csv, selection policy, coverage curve, proposed references |
| Same-pocket consensus | structure_consensus | Chosen ligand instance/target chain, protein/pocket alignment, admission and prepared templates |
| Consensus proposal/design | consensus_recommend, optional adopt/design | Actual contact evidence, supported template states, unresolved review requirements |
| Budget screening | consensus_budget | Connected library, retrieval and output molecule budgets, scope and actual completion |
| Additional delivery | budget_page | Intended ranking branch and latest delivered rank interval |

Initial same-site contact overlap in the diversity stage is not the final aligned
pocket admission. Chemical diversity is not 3D pose diversity or activity evidence.
Polymer ligands retain their separate chemical-preparation limitations.

Budget screening retrieves from the whole configured USRCAT FAISS library index,
expands stored conformers for retrieved molecules, and performs the existing
same-pose scoring and contact-first template ranking/RRF. It does not score every
library molecule exhaustively. Threshold-based `guided` is a different branch.
The current target's ANN recall and biological enrichment remain validation tasks.

## English prompts

Use the conversation containing your completed tasks. A new conversation does not
implicitly inherit another session's task IDs; attach its existing owned plan using
the existing import control if needed.

Start by locating work already done:

> Inspect my structure-to-library workflow. Identify the PDB collection, ligand diversity, aligned consensus, recommendation and completed screening tasks. Show their source relationships. Do not start new computation.

For a new target, substitute the actual accession and reference PDB:

> For human MDM2, UniProt Q00987, use reference PDB 6Q9L to collect experimental ligand-bound structures and compare ligand diversity. Report the quality and site exclusions, similarity matrix, coverage and proposed diverse references. Do not search the molecule library yet.

Review existing selection:

> Explain the ligand differences in task TASK_ID using its actual similarity matrix and coverage results. Which references add chemical diversity, and which are redundant? Distinguish chemical similarity from aligned 3D differences.

Build a consensus with actual identifiers confirmed in the report:

> Build the aligned same-pocket consensus from diversity task TASK_ID using reference ligand 6Q9L:HTZ:A:201 and target chain A. Report admission failures and the prepared template set.

> Generate the consensus recommendation from task CONSENSUS_TASK_ID. Explain the contact evidence, selected templates and unresolved scientific limitations.

Then explicitly request the budgeted search:

> Use recommendation task RECOMMENDATION_TASK_ID for budgeted screening of my connected library. Retrieve 1000000 unique molecules, rank contacts first with shape as a tie-breaker, fuse template ranks and export the first 100000 unique molecules with original names and original and posed MOL2 files. Use my configured workers and chunk size. Report any shortfall and the actual search scope.

After completion:

> Read the completed screening task SEARCH_TASK_ID. Show the exported molecule count, ranking scope, file locations and chemistry review flags. Do not rerun the search.

> Export the next 100000 unique molecules after the last delivered page of this search. Reuse its saved ranking without rescoring, and report the rank interval.

View the reference structures:

> Open the aligned consensus structures from task CONSENSUS_TASK_ID in PyMOL and show the ligand and nearby whole residues. Explain which references each object represents.

## Execution and acceptance

Update the existing branch and restart Chat with the same workspace/runtime/provider
configuration as described in E072. Stages remain asynchronous scientific tasks.
After a task completes, a follow-up such as `Continue to the next authorized stage`
lets the model inspect actual dependencies and choose the next executor. This change
does not install an unattended whole-chain scheduler. It does not infer authorization
to launch an expensive search from a request to inspect results.

The tool offers no next action for failed child reports or unresolved reference
readiness. Choose from actual reference options; do not invent a chain or ligand ID.
When multiple targets/branches match, select the intended task. Missing source tasks
retain their source-run IDs instead of being assigned to an unrelated branch.

Local regression covers source lineage, similarity-file reads, reference ambiguity,
failed reports, session separation and path containment. Live PDB retrieval, model
choice quality and the workstation screening pipeline still require acceptance on
the configured system. No full-library search was launched for these local tests.
