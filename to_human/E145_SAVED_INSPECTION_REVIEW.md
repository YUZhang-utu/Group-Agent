# Saved PDB inspection review

The E144 inspection may finish successfully even when a subsequent model-driven
report summary fails JSON validation. Do not rerun the scientific task.

E145 adds `pdb_review` and a narrow natural-language route for an explicit single
inspection task ID. It reads the owned saved report, verifies report and structure
digests, renders deposited facts and records display coverage without a provider
call. It performs no remote retrieval, new calculation or receptor adoption.

## Workstation acceptance

Pull and restart Chat using the same storage. In the owning target project, send:

```text
Read the saved inspection report for completed task TASK_ID.
Summarize its protein identity, deposited mutation evidence, chain coverage,
biological assembly annotations, ligands and ions, database mapping and saved paths.
Record a concise project note linked to this task.
Do not repeat the inspection or launch any calculation.
```

Replace TASK_ID with the completed task's actual ID. A successful trace records
`answer_mode: local_pdb_report_review_no_model_call`. No provider is required for
this route. The response shows note status and its saved ID; repeating the same
request or opening another chat in the same project reuses the note. Without an
explicit note request, only a read audit is written. Different projects cannot
read the task. Changed or missing report/structure files stop the review.

The display includes at most 20 rows per section and marks shortened fields.
Full sequences, residue lists and assembly operators remain in the saved report;
the summary does not claim to show every field. Deposited mutation/mapping entries
are not independent isoform validation. Missing metadata is not evidence of absence.

The broader report-analysis planner, mixed compute requests, literature synthesis
and generic blocked-plan clarification are unchanged. This is a local PDB-report
review path, not a claim that all model-output failure modes have been solved.
