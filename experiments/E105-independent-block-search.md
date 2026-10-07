# E105: independent leading-block search comparison

Protocol authorized 2026-10-07; implementation and local fixtures precede workstation execution.

Compare E094/E095/E096 x ChemPLP/EquiScore x Top-5/Top-10 blocks (12 arms), using the existing 7NA2_A ranking context and identical confirmed MDM2 ligand query. Block ranking remains the mean of the best ten distinct molecules. Do not mix scoring arms or average their score units.

Process blocks in ascending saved rank. Search every conformer of each started block with the existing contact-first refinement funnel. At each completed-block boundary, stop if cumulative qualified conformers reach 100000; rank the searched prefix with existing per-template conformer RRF and export at most 100000. Otherwise continue up to the selected block limit and report shortfall. Multiple conformers per molecule are allowed. Top-5 and Top-10 can be identical if both stop within the first five blocks.

No docking or rescoring is authorized for these arms. Preserve source MOL2 chemistry and record source identity and search transforms. Exported MOL2 files contain original coordinates, not newly docked poses. Record actual visited blocks, searched/ranked/exported counts, unique molecules, per-block timing, shortfall and hashes. Resumed timing is not a fresh speed benchmark. Candidate counts and scores are not experimental enrichment or affinity.

Preserve the completed E095 union search/docking run PROMPT-874365d177704bf3 as a separate historical result (user receipt: 344476 searched, 100000 docked, 72859 molecules, zero failed jobs). Do not reuse its truncated union export as either independent arm or overwrite its sealed outputs.

Validation: fixture tests for ordered stopping, undersupply, independent arm dispatch, no docking dependency, owned inputs, idempotent requests, receipt integrity and MOL2 export. Real-library execution is pending workstation launch; no speed or enrichment claim yet.
