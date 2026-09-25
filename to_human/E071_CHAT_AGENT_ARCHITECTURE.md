# E071: composable tools and research-agent gaps

User requirement, 2026-09-25: understand natural requests, operate PyMOL, inspect
existing results, consult research literature, and produce grounded explanations.
Adding a phrase-specific branch for each user request is not the desired design.

## Current repair

- Expose the existing pocket_view operation to the PyMOL planner as a composable
  adapter: objects is a space-separated list of catalog IDs; radius is 1..12 A.
- Preserve catalog ligand identity and whole residues, enable selected objects,
  label residues, export the residue list, and keep program checkpoint/rollback.
- Pass bounded recent conversation into the PyMOL planner; previously only the
  outer router received that context.
- Include the final underlying error in the top-level failure message.
- Accept pymol_agent with or without the leading slash. No phrase-specific
  interpretation bypass is introduced; the model selects and composes tools.
- Simple pocket programs no longer load unrelated aromatic metadata.

Workstation example after updating and restarting Chat plus bridge PyMOL:

```text
/pymol_agent Show residues within 5 angstroms of the ligand in v001. Label the residues and keep the ligand visible.
```

Existing direct operation for diagnosis:

```text
/view {"operation":"pocket_view","objects":["v001"],"radius":5}
```

The locally tested planner is a fixture, not the user's live provider. The user's
generic two-attempt failure message does not identify the original failing call.
Local API-contract tests do not prove arbitrary natural-language comprehension.

## Required general architecture: not yet implemented

One conversation coordinator should maintain the user's objective and evidence,
then repeatedly choose a tool, inspect its receipt and decide the next action or
answer. The tool registry should expose capabilities, not phrase-specific routes:

1. Inspect current scene, selections and measured geometry.
2. Compose display operations and validated scientific analyses.
3. Read owned task reports, tables, provenance and confidence summaries on demand.
4. Search literature and retrieve source text with bibliographic provenance.
5. Synthesize an answer distinguishing measured results, source claims and hypotheses.

Read-only tool calls and successful execution receipts must feed back into the
same model conversation. Use bounded steps and time, explicit pending states,
resumable execution and action-specific permissions. Do not blindly retry a job
that may already be running. A successful API call is not proof of scientific
validity or visual correctness. Unsupported analyses must explain missing tools
or inputs rather than fabricate a result.

Acceptance should include paraphrases and multi-turn requests: show the pocket,
hide unrelated chains while retaining the ligand, explain computed contacts,
compare two completed results, find supporting/contradicting papers, cite those
papers and identify limitations. Record real-provider tool traces and reference
outputs. Literature search and a general research tool loop are absent from the
current PyMOL adapter; bundled skill text does not supply these capabilities.
