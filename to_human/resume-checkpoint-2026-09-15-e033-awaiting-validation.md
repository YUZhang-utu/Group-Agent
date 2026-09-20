# AIDD progress checkpoint: 2026-09-15

Historical record, superseded by later E033/E034/E035 results. At this date the
user was still waiting for registration/precomputation; E033 had not been validated.
Local code tests were not evidence of real-library acceptance.

Project: `D:/agent/projects/aidd_macrocycle_agent`.
Remote: `https://github.com/YUZhang-utu/Group-Agent`, branch main.
Latest pushed commit at the time: `bbd7675`.

## Completed at this checkpoint

- `bab0643`: retry random compact database-ID collisions.
- `14e156d`: read-only conformer-conflict diagnosis and explicit error location.
- Workstation `split_0119.mol2` records131568and131625(zero-based) both named
  `c--L-dA-Lnme-Wnme-dL-VNMe-c_conf1`, with116atoms/118bonds, identical topology hashes
  but different coordinates. This was not the earlier random ID collision;
  stereochemical equivalence was unconfirmed.
- `65317f6`: `--preserve-index-conflicts` keeps distinct same-index records with
  separate conformer IDs and negative internal indices, preserving original names,
  source positions and audit information without changing source files.
- `1e94c8c`: E033 protocol; `80acb64`: acceptance/evaluation; `bbd7675`: sampling notes.
-147local tests passed/2dependency skips; Bash syntax passed. Real FAISS/RDKit and
  workstation validation remained pending.

## Historical next execution

Repository: `/mnt/medchem_taltio/wrk/yu_agent/Group-Agent`.
Precompute: `/mnt/local/hand/yuzhang/aidd/library-precompute-20260911`.
Registration, features, training and merging must finish with `COMPLETE.json`.
Wait for an active process; otherwise resume the original command in
`E032_CONFORMER_CONFLICT.md` with `--preserve-index-conflicts`.

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
bash scripts/run_e033_library_acceptance.sh
```

The runner checks completion and acquires the existing batch lock without changing
the source library. Results use a sibling timestamped `e033-library-acceptance-*`
directory. See `E033_LIBRARY_ACCEPTANCE.md`.

Pending evidence at the time: E032completion/integrity; E033report.md/report.json
(and metrics/failure logs if needed); actual latency, global exact-USRCAT
Top100/1000recall and conformer/molecule retention by budget. Do not replace absent
measurements with estimates. Coarse descriptor calibration is distinct from
activity enrichment or refined pose quality. No new experiment was launched merely
to save this checkpoint. Use current research-state.yaml and newer result records
when resuming today; do not treat this historical waiting state as current.
