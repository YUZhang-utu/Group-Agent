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
