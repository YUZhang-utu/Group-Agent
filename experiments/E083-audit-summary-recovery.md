# E083: audit summary disk failure recovery

The workstation failed in the ring-size GROUP BY after conformer inserts were
committed. Hypothesis: summary sorting can exhaust SQLite temporary storage even
when the database has space. The specific filesystem cause remains unconfirmed.

Protocol before checks: replace global JSON histogram sorts with streaming
counters, preserve summary semantics against reference SQL, and provide an
explicit summary-only recovery mode. Recovery must reject incomplete source
manifests, changed source hashes and database/source count mismatches. Preserve
the original report, issues, scientific code hashes and database bytes. Test a
simulated post-ingestion failure, refusal cases and retry behavior. No chemistry
recomputation or clustering acceptance claim is made by this engineering repair.
