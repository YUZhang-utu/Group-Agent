# E116: shared KRAS G12D / G12V Switch-II screening

User-confirmed scope: noncovalent Switch-II screening on 24 CPU cores / 64 GB RAM.
The runner uses 20 processes total, one BLAS thread each, with one library decode
per conformer for both branches. All catalog conformers are scheduled in full
mode. No top-block exclusion, ranking cap or docking is substituted for screening.

## Hard screening policy

All required conditions must hold in the SAME pose for its paired native receptor:

| Condition | Both branches | G12D-specific | G12V-specific |
|---|---|---|---|
| His95 | N/O acceptor-compatible polar proximity, <=3.5 A | Required | Required |
| Asp69 and Glu62 | N/O donor-compatible polar proximity, <=3.5 A, both regions | Required | Required |
| Met72 | Hydrophobic/aromatic carbon to sidechain carbon, <=4.5 A | Required | Required |
| Arg68 | Polar acceptor/guanidinium or aliphatic carbon contact | Polar mode required | Either reference-compatible mode |
| Residue 12 | Evaluated independently | N/O donor to Asp12 oxygen, <=3.5 A | Val12 steric compatibility; no invented obligatory H bond |
| Native feature geometry | Directional one-to-one score >=0.5 for required polar modes; >=3 distinct matched features/residues | Asp12 and Arg68 polar modes included | No Asp12 feature imposed |
| Whole shape/color | Gaussian point-cloud Tanimoto >=0.50 / >=0.30 | Same pose | Same pose |
| Sterics/coverage | Minimum heavy-atom distance >=2 A; native contact envelope >=0.65; each of three spatial regions >=0.50 | Native receptor | Native receptor |

These are strict, exploratory operational gates. They are not proven necessary
conditions for all KRAS inhibitors, validated H-bond angles, or affinity estimates.
No interaction energy, His tautomer enumeration or water bridges are modeled.
Organic fluorine cannot fulfill mandatory N/O polar conditions. The whole-color
baseline remains unchanged; E115 experimental scoring is not adopted. Output
near misses are diagnostic only and never enter the passing-candidate CSV.

## References and limitations

G12D: 7RPZ/6IC (MRTX1133). G12V: 12FF/A1DB0 (macrocyclic compound 20).
Each uses its own GDP-bound receptor plus same-chain GDP/Mg. Water and other
crystallization components are excluded. Blank/A alternate atom locations are
selected with highest occupancy per identity. Both native constructs include
C51S/C80L/C118S engineering backgrounds; G12D also has an observed N-terminal
identity/numbering difference. No sidechains were silently repaired. This is
experimental-construct screening, not proof of canonical-mutant selectivity.
One reference chemotype/state per allele limits scaffold/state coverage.

Reference package: configs/targets/kras_sip_e116. The package records coordinate,
CCD and KRAS4B sequence hashes, atom/feature metadata and all screening thresholds.
Native self-recovery passes for both references. Shared versus independent
evaluation, chemical negatives, two-process execution, deduplication, resume and
corruption checks are covered by local tests (31 focused checks passed, including
eight E116 tests and actual two-process fixture execution). Reference seals also
pass Windows/Linux newline/path portability checks. Workstation timing/pass rates are
pending; synthetic tests do not establish library enrichment.

## Run today

Use the existing activated aidd_workstation environment and repository root.

```bash
git pull --ff-only origin feature/structure-guided-chat
bash scripts/run_kras_weekend.sh pilot
cat "$HOME/kras-e116/pilot/report.json"
```

The pilot examines 1000 spread conformers against BOTH targets, with a one-hour
cooperative cap (normally intended to be much shorter, actual time unmeasured).
It reads the library location from $HOME/mdm2-e109/prepared/report.json. If that
receipt is elsewhere, pass --library-receipt /actual/prepared/report.json; or pass
--batch /actual/library/root containing artifacts/catalog.json and
chemical/catalog.json. Pass the same location for pilot and full runs.

Inspect full-pass fractions, per-stage rejection and the linear dual-target
runtime estimate. Zero pilot hits is not proof of no library binders. High pass
fractions require reviewing whether the conditions are selective before full
launch; do not simply tune thresholds to achieve an arbitrary rejection target.

After reviewing the pilot, launch the full search in a persistent terminal:

```bash
tmux new -s kras-weekend
bash scripts/run_kras_weekend.sh full
```

Detach with Ctrl-B then D. This starts both targets together, not 20 processes
per target. The 60-hour cooperative budget stops new chunks and finishes in-flight
chunks. It does not guarantee completion by Monday. Re-run the same command to
resume matching sealed chunks. Do not update code, reference files or library
artifacts while a run is active. Input/code mismatch refuses reuse.

Optional: --hours 48 changes the invocation budget; --workers 16 reduces resource
usage. Chunk size and sampling settings must remain unchanged on resume. Keep
the machine awake and ensure output storage is available. No new environment,
GPU package or API key is required for the scanner.

## Results

- progress.json: current per-target conformer coverage and stage counters.
- report.json: complete/partial status, explicit full_library_complete flag,
  unique passing molecules and runtime estimate.
- KRAS_G12D_candidates.csv / KRAS_G12V_candidates.csv: passing molecules only,
  deduplicated by stable molecule ID, with best passing pose and reference ID.
- candidates.sqlite: rebuildable distinct-molecule candidate table.
- chunks/*.json.gz plus receipts: all passing conformers and bounded near misses.
  Sealed chunks are the recovery source; keep them for later export and audit.

Counts: chemical is shared coarse eligibility; contact includes direct native
region and one-to-one feature gates; shape and color are independently counted
after contact; geometry additionally requires shape; full requires all on one pose.
Ranks use shape+color only after hard acceptance, not as a substitute for it.
CSV output is an ID/pose handoff, not prepared docking ligands or binding results.

For an initial return, include pilot/full reports, progress, controls, both CSVs
and logs. Preserve all chunks on the workstation; they may be too large to upload.
Candidate SDF/MOL2 export and docking preparation are the next separate step.

Sources: https://www.rcsb.org/structure/7RPZ,
https://www.rcsb.org/structure/12FF,
https://rest.uniprot.org/uniprotkb/P01116-2.fasta.
