# AIDD two workstation jobs completed - 2026-10-01

Evidence: user-pasted completion summaries and E094 validator receipt. No direct workstation filesystem access or independent remote rehash performed.

## E094 result
Output: /mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/final-blocks-e094
Structural and property gates passed. Total 24,663,736; regular 24,530,851; special 132,885 (0.538787%). Parent work blocks 385 including special pool; property blocks 384 regular. Minimum 5,393; none below 5,000; no logical capacity ceiling.
The implementation refines each regular parent independently without merging parents: 384 regular parents and 384 property children imply no added splits. Profiles were used for refinement, but this does not demonstrate that property similarity is scientifically calibrated. Inspect properties/blocks.csv split_evidence and profiles/report.json before considering threshold changes. Do not force extra splits or rerun the full build on aggregate counts alone.
Average regular block size is approximately 63,882.4. Inspect full size distribution, largest blocks, mixed-class composition and split evidence. Next scientific evaluation: stratified 3D retrieval, PLANTS/MDM2 and user N-E rescoring; preserve exhaustive special-pool search and do not reject unsampled blocks. Chemistry review population 1,150,072 remains separate and deferred.

## Top-million conformer export
Output from launch handoff: /mnt/local/hand/yuzhang/aidd/exports/mdm2-11templates-top1m-conformers
User summary: status complete; input/scored 3,854,146; unscored 0; exported 1,000,000; distinct molecules 460,332.
Expected default part size: 10,000; top-00001.mol2 through top-00100.mol2 with same-name CSV and JSON receipts. Actual part inventory is authoritative in report.json parts_manifest. selected-conformers.csv is the combined ranking/identity/file-location manifest; selection.sqlite retains ranking state.
Coordinates are original library conformers, not docking poses. Ranking uses saved scored ANN candidates, not exhaustive full-library affinity scoring. Multiple conformers per molecule are intentional.
Export is already done: transfer MOL2 parts plus matching CSV/JSON, selected-conformers.csv and report.json using the workstation file transfer method. No selector rerun or merging is required. selection.sqlite need not be copied solely to consume MOL2 structures.

## Resume
The running-job checkpoint is superseded for status only. Both jobs are user-reported complete. Preserve outputs and provenance. Next obtain properties/blocks.csv, backbone/report.json and profiles/report.json for diagnostic review; locate and transfer the completed top-million export. Do not change MOLIQ or trigger new heavy calculations automatically.
