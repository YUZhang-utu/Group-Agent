# E032b — Auditable MOL2 kekulization fallback

## Status

Confirmatory protocol locked before implementation and test execution.

## Problem

Some new macrocycle MOL2 records carry explicit Tripos aromatic atom and bond
types that RDKit can preserve but cannot convert to one consistent Kekulé
single/double-bond assignment. `MolFromMol2Block(..., sanitize=True)` therefore
returns `None`, and one multiprocessing worker aborts the complete E032 shard.
The same strict parse is repeated by both artifact-v1 and chemical-companion
generation.

## Hypothesis

If strict parsing fails, parsing without sanitization and then running every
RDKit sanitization operation except `SANITIZE_KEKULIZE` will admit precisely the
records whose only unsupported operation is Kekulé assignment. Their original
aromatic graph and coordinates should remain usable for BaseFeatures, PMI,
USRCAT, chemical topology, feature directions, and terminal-torsion discovery.

## Intervention

- Add one shared MOL2-to-RDKit loader.
- Keep strict full sanitization as the primary path.
- Permit a fallback only when all sanitization operations other than
  `SANITIZE_KEKULIZE` succeed.
- Reject records that fail the reduced sanitization set.
- Use the shared loader in both conformer-artifact and chemical-companion
  workers.
- Record strict/fallback counts in both shard manifests.

## Prediction

- A synthetic five-member Tripos aromatic ring that reproduces the RDKit error
  is accepted in `aromatic_no_kekulize` mode.
- Ordinary MOL2 remains in `strict` mode.
- Structurally invalid MOL2 remains rejected.
- The supplied representative failing macrocycle yields features, PMI, and a
  60-value USRCAT vector under the fallback.
- Existing tests remain green.

## Decision rule

Accept the change only if targeted parser/artifact/companion tests pass and the
representative real record passes every non-kekulization sanitization operation
and downstream descriptor smoke checks. Linux pilot on `N5_0.mol2` remains the
required workstation confirmation before resuming the full frozen batch.

## Confirmatory result — 2026-09-14

The Linux pilot passed under Python 3.11.16, RDKit 2026.03.5 and FAISS 1.14.3.
The driver preserved the failed partial, rebuilt `N5_0.mol2`, and completed all
5,744 inserted conformers with zero duplicates in 6.7 seconds. Artifact-v1 and
chemical-companion manifests independently reported identical sanitization
counts: 4,485 strict and 1,259 `aromatic_no_kekulize`. Both sums equal 5,744;
there was no record loss or cross-stage disagreement. E032b is confirmed for
this shard. The frozen full-library E032 run is now the next action.
