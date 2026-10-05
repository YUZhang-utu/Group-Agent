# Published rescoring models for manual review

Checked 2026-10-05 against publisher pages and author repositories. Selection
prioritizes rescoring saved PLANTS complexes and early enrichment, not generating
new poses. User declined GNINA as the desired public model. No installation,
inference or automatic model selection is performed. This replaces the earlier
GNINA-first proposal, while preserving the user-owned N-E interface.

## Primary shortlist

| Model | Publication | Why inspect it | Main review question |
|---|---|---|---|
| EquiScore | Nature Machine Intelligence 6, 688-700 (2024) | External-pose screening with explicit cross-docking-engine evaluation; author code provides a screening entry point and weights | Does its input representation retain our noncanonical macrocycle chemistry? |
| SCORCH2 | Advanced Science 12, e08318 (2025), online August 20 | Interaction-feature consensus targeting enrichment, with external docking poses | Inspect train/test overlap, deduplicated variants and dependence on docking source |
| GenScore | Chemical Science (2023), D3SC02044D | Extends RTMScore to balance scoring, docking, ranking and screening | Select and freeze the appropriate scoring variant rather than assuming all checkpoints are interchangeable |
| RTMScore | Journal of Medicinal Chemistry (2022), 2c00991 | Residue-atom distance statistical potential with pose and screening benchmarks | Inspect pocket construction and large-ligand coverage; this is an older baseline |
| DeepRLI | Digital Discovery (2025), D4DD00403E | Distinct scoring, docking and screening readouts | Use the screening readout for enrichment; keep affinity/pose outputs separate |

Primary-source links, including full titles:

1. EquiScore: *Generic protein-ligand interaction scoring by integrating physical
   prior knowledge and data augmentation modelling*.
   [Paper](https://www.nature.com/articles/s42256-024-00849-z) |
   [Code](https://github.com/Intelligent-Drug-Discovery-Lab/EquiScore).
2. SCORCH2: *A Generalized Heterogeneous Consensus Model for High-Enrichment
   Interaction-Based Virtual Screening*.
   [Paper](https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202508318) |
   [Code](https://github.com/LinCompbio/SCORCH2) |
   [Weights/data](https://zenodo.org/records/14994007).
3. GenScore: *A generalized protein-ligand scoring framework with balanced
   scoring, docking, ranking and screening powers*.
   [Paper](https://pubs.rsc.org/en/content/articlehtml/2023/sc/d3sc02044d) |
   [Code](https://github.com/sc8668/GenScore).
4. RTMScore: *Boosting Protein-Ligand Binding Pose Prediction and Virtual
   Screening Based on Residue-Atom Distance Likelihood Potential and Graph
   Transformer*.
   [Paper](https://doi.org/10.1021/acs.jmedchem.2c00991) |
   [Code](https://github.com/sc8668/RTMScore).
   The newer [Nature Protocols workflow](https://www.nature.com/articles/s41596-026-01389-z)
   was published June 24, 2026; it integrates RTMScore, not a new 2026 RTMScore model.
5. DeepRLI: *A multi-objective framework for universal protein-ligand interaction
   prediction*.
   [Paper](https://pubs.rsc.org/en/content/articlehtml/2025/dd/d4dd00403e) |
   [Code](https://github.com/fairydance/DeepRLI).

## Adjacent recent candidates

- DeepPpIScore: *Harnessing deep statistical potential for biophysical scoring
  of protein-peptide interactions*. Acta Pharmacologica Sinica 47, 518-532
  (2026 issue; online October 1, 2025).
  [Paper](https://www.nature.com/articles/s41401-025-01659-8) |
  [Author-linked code](https://github.com/zjujdj/DeepPpIScore).
  Peptide-focused geometry/statistical scoring makes it worth inspecting for our
  chemical class. Direct support for noncanonical residues, N-methylation and
  ring closure is unverified. The code link was obtained from the paper; its
  contents were not accessible in this browsing pass.
- E-CloudBind: *An electron-density point-cloud framework for robust protein-ligand
  interaction prediction*. Nature Communications 17, 7424, June 11, 2026.
  [Paper](https://www.nature.com/articles/s41467-026-74196-5) |
  [Code](https://github.com/Liuyujian0408/DPI).
  Recent structure-based affinity regression with electron-density preprocessing.
  Treat as an affinity-side candidate, not established evidence of superior pose
  rescoring or EF1 on PLANTS macrocycles. Audit preparation cost and pose retention.

## Interpretation and proposed review order

Read EquiScore and SCORCH2 first for the current enrichment-oriented task; GenScore
is the next broad scoring candidate. Read DeepPpIScore alongside them for peptide
relevance. These are project-fit judgments, not a common-benchmark accuracy ranking.
Digital Discovery is included for specialist relevance, without claiming its
journal standing equals Nature Machine Intelligence or Chemical Science.

SCORCH2 Table 1 provides a useful same-pose illustration: reported EF1 for Glide
is 12.4722, Glide plus EquiScore 16.832, and Glide plus SCORCH2 19.264. Other docking
sources change the relative order. This is author-reported benchmark evidence,
not an expected improvement for our target. The paper discusses training overlap
and fully deduplicated variants; inspect those before accepting generalization.

None of this source review establishes a strongest model for our particular
noncanonical macrocycle library. Audit atom typing, ring/stereo retention, size
dependence, training overlap and early enrichment on the same labeled test panel.
Without activity labels, rank agreement remains exploratory. Preserve model
native direction/units: higher-is-better outputs cannot use ChemPLP's ascending
sort. After rescoring, recompute each molecule's best conformer and each block's
Top-10; do not merely rescore the ten ChemPLP winners and call it a full reranking.
