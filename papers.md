# Papers and software references

## Reduce / hydrogen placement

- Word JM, Lovell SC, Richardson JS, Richardson DC (1999), *Asparagine and
  glutamine: using hydrogen atom contacts in the choice of side-chain amide
  orientation*, Journal of Molecular Biology 285:1735--1747. Richardson Lab
  documentation describes Reduce `-build` optimization of adjustable hydrogen
  groups and Asn/Gln/His side-chain orientations. Project relevance: query-only
  pocket hydrogen/flip validation before reporting hydrogen-bond angles; never
  a library recall filter. https://github.com/rlabduke/reduce

## E033 engineering references — 2026-09-15
FAISS official wiki: https://github.com/facebookresearch/faiss/wiki/How-to-make-Faiss-run-faster and https://github.com/facebookresearch/faiss/wiki/FAQ . Used to ground search effort, candidate budget and timing distinctions; our .95 gate and panel are project protocol choices, not source guarantees. No biological filtering rate is inferred from these sources.

## E038 engineering API references — 2026-09-19
Official API documentation (not research papers): OpenAI Chat Completions https://developers.openai.com/api/reference/resources/chat ; UniProt REST query help https://www.uniprot.org/help/api_queries ; RCSB Data API https://data.rcsb.org/ ; AlphaFold3 inputs https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md . Used for provider-neutral JSON chat, public sequence/structure evidence, and basic AF3 dialect-v1 protein/CCD inputs compatible with existing local installations. No new biological findings from these API references.


## E046 implementation references - 2026-09-21

CuPy official installation: https://docs.cupy.dev/en/stable/install.html ; CUDA
runtime component wheels and version-specific distributions. CuPy official
performance guidance: https://docs.cupy.dev/en/stable/user_guide/performance.html ;
asynchronous execution, synchronization and first-use costs. Engineering sources,
not evidence of biological screening quality or measured RTX 5090 acceleration.

## E051 sources

- RCSB Search/Data API: https://search.rcsb.org/ and https://data.rcsb.org/ .
  Target-associated PDB discovery and entry metadata; metadata alone do not
  establish prepared receptor quality or ligand binding relevance.
- Europe PMC REST API: https://europepmc.org/RestfulWebService . Search/core
  abstracts are provenance-bearing leads, not automatic experimental validation.
- Pharmit help: https://pharmit.csb.pitt.edu/help.html . Pharmacophore feature
  selection and receptor exclusion inform workflow design. E051 does not claim
  Pharmit algorithm equivalence or adopt its thresholds as a calibration.

## E052 source evidence

- UniProt Q00987: https://rest.uniprot.org/uniprotkb/Q00987.json . Verified human
  MDM2 sequence and accession, retained locally with source payload checksum.
- RCSB entry: https://www.rcsb.org/structure/5C5A . Nutlin-3a complex supplies the
  project reference pocket; original coordinates supersede missing export metadata.
- RCSB entries https://www.rcsb.org/structure/7NUS ,
  https://www.rcsb.org/structure/8GCG and https://www.rcsb.org/structure/9CDZ .
  Deposited polymer entities establish why peptide ligands require a separate path.
- Official Data API: https://data.rcsb.org/ . Cached Search/GraphQL/CCD responses
  and raw coordinate files under the project analysis directory record actual counts.

## E082 MDM2 ensemble docking and PLANTS

Primary-source notes: literature/e082-mdm2-ensemble-sources.md. Relevant evidence:
Bista et al., Structure (2013), doi:10.1016/j.str.2013.09.006; MDM2 ensemble receptor
models, JACS (2007), doi:10.1021/ja073687x; PLANTS scoring, JCIM (2009),
doi:10.1021/ci800298z. These motivate an ensemble/cross-docking protocol, not a
claim of measured performance for the user's macrocycle library or N-E workflow.

## E084 pocket-state evidence

POVME 3.0: Software for Mapping Binding Pocket Flexibility:
https://pmc.ncbi.nlm.nih.gov/articles/PMC5751414/ . Pocket shape/chemical comparison
motivates the workflow; our bounded grid implementation is not a POVME replica.
MDM2 transient-state experimental evidence:
https://pmc.ncbi.nlm.nih.gov/articles/PMC4104591/ . Rare states warrant review;
PDB deposition counts are not equilibrium occupancy estimates.

## 2026-09-29 E088 descriptor reference

RDKit official documentation:
https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-3d-descriptors
Accessed for distinction between chemical properties and coordinate-dependent
shape descriptors. E088 local VDW-axis extents and branching are explicitly
implemented proxies, not a claim of exact volume or docking-energy validation.

## E092 clustering interpretation references (2026-09-29)

U. von Luxburg, Clustering Stability: An Overview:
https://arxiv.org/abs/1007.1075
Stability has limitations as a universal selector of true cluster count.
U. von Luxburg et al., Clustering: Science or Art?, PMLR27 (2012):
https://proceedings.mlr.press/v27/luxburg12a.html
Cluster evaluation depends on task context. Our proposed PDB pocket evidence
gate is a project design, not an established physical continuum test.
