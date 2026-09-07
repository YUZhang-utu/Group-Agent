# E026 Real 8BJU/1X8B multi-cocrystal validation protocol

Status: protocol locked before workstation execution

## Fixed inputs

- Query 1: 8BJU, ligand QT9, the already accepted real detailed refinement.
- Query 2: 1X8B, ligand CCD 824, auth chain A, residue 901.
- Biological site: `WEE1-ATP-site` for both queries.
- Receptors remain distinct: `8BJU-prepared-v1` and `1X8B-prepared-v1`.
- Library: `LIB-AFA68EE6888C`, 100,000 molecules and 299,999 conformers.
- Primary aggregation objective: `atomcentered_anchored_joint`.

RCSB describes 1X8B as the human WEE1A catalytic domain complexed with the
active-site inhibitor PD0407824 at 1.81-A resolution. The deposited ligand is
CCD 824 at auth chain A/residue 901. These identifiers must be verified from
the downloaded mmCIF before scoring.

## Hypothesis

An independent 1X8B/824 search will yield a valid detailed Gaussian result and
E025 will form the site-scoped union with 8BJU/QT9 without losing unique-query
molecules. Shared molecules will receive two-query support and two receptor-
specific docking tasks; unique molecules will retain one supported task.

## Locked workflow

1. Download official 1X8B mmCIF and CCD 824, retaining file hashes.
2. Verify exactly one 824 instance resolves at auth A:901.
3. Extract distance/SASA-supported polar anchors and compile the query plan.
4. Reuse the existing library artifact, FAISS and pharmacophore indexes.
5. Run the complete tiered retrieval validator with FAISS baseline preservation.
6. Prepare the Gaussian query and run resumable PCA coarse plus pair refinement.
7. Build a two-query E025 plan with the accepted QT9 and new 824 refined results.
8. Aggregate by molecule and emit receptor-specific docking tasks.

## Confirmatory acceptance criteria

- Downloaded files parse and 824 resolves uniquely to A:901.
- Pharmacophore validation reports `accepted=true`; every FAISS baseline ID is
  retained. Zero supported anchor pairs is permitted but must be explicit.
- Gaussian run reports `status=complete`, finite objective scores, valid hashes,
  and no missing selected IDs.
- E025 reports exactly two queries and one site.
- The summary key set equals the union of per-query molecule key sets.
- At least one molecule has support from each query; zero shared molecules is a
  reportable scientific result, not an orchestration failure.
- Every shared admitted molecule has both receptor-specific docking tasks.
- Every unique admitted molecule has its supported receptor task.
- Raw cross-query scores remain unaveraged and all E025 invariants pass.

## Classification boundary

This experiment validates multi-query orchestration and measures overlap for
two real WEE1 co-crystals. It does not establish enrichment, binding, docking
accuracy, or superiority of either query.

