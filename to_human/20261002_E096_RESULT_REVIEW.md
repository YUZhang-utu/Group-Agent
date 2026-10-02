# E096 workstation receipt review - 2026-10-02

## Evidence received

User supplied completion reports, not local copies of the large output databases.
Profiles complete: 24,530,851 / 24,530,851, coverage 1, width 317.
Regular blocks: 1,189 from 384 original E094 parents; 805 accepted splits.
Minimum 5,015; maximum 94,952; zero blocks below 5,000. Maximum children 4,
minimum original-parent fraction 0.20. Special pool unchanged: 132,885, about
0.538787% of 24,663,736 admitted conformers. The separate chemistry-review
population remains deferred and is not part of this percentage.
Membership receipt: passed_against_complete_cached_profiles.

Feature activity: chemistry 19/33; local sterics 33/36; backbone geometry 22/24;
typed spatial 162/224. Reported aromatic/hydrophobic counts are nonzero.
Positive/negative ionizable residue-contained feature counts are zero.
Check grouping is source molecule ID CRC32 modulo five.
Extraction invocation: 37,917.5816 s (10.533 h); partitioning invocation:
1,285.2944 s (21.422 min). These are reported stage times, not independently
observed whole-job wall-clock duration.

User-reported hashes:
- profiles.sqlite: 63e6f98959cfb82c6aca842087b7984b8d6fe00806402f5cffe9d5d72bab2f8a
- property_blocks.sqlite: 59595021088e8d1df3c918780214bbf3f6667f4dc79b5212f561688ea6999b69
- blocks.csv: eb4f583fae2d5b7c757fb89c5a2cd89e2ea5d39ab7e930f4f73ca2baf2d950a1
- transform.json: cca054c50ee8266baf3581f979a84baa5487db0c75dd39cc6ecc6b6c312318c2
- split-decisions.json: 3ded272a1e76b0bf1b09b1e09097360609c6d38c6df29a19e0ccd8de47d4f2cf

## Interpretation and decision

Reported completeness and fragmentation gates pass. Retain E096 as the current
offline joint-block candidate; do not tighten splitting or rerun source extraction.
E095 had 1,188 blocks (minimum 5,011; maximum 94,822). Similar counts are unsurprising
under the same parent/floor/four-child policy; they do not prove identical memberships
or lack of spatial influence. The only reported rejection reason is the population
floor; this also means count alone is a weak measure of feature quality.

Nonzero active geometry/spatial channels establish variation, not screening quality.
Nominal group budgets need checking against weighted variance after scale floors.
Splitting-axis coefficients can show participation but are not causal importance.
Source molecule checks remain internal/adaptively reused, not unbiased recall tests.

## Ionizable-channel investigation

The code activates all seven RDKit feature families and maps residue-contained
features into spatial descriptors. Four synthetic cyclic peptide controls carrying
neutral/protonated amines and acid/carboxylate groups passed the local descriptor
mapping tests. However, local RDKit 2026.03.6 has feature-definition hash
849316d651bf499859058aa82680a59bc3fdb610c35e3163518448bdfffdf782,
whereas the reported run used 2026.03.5 and
766f9790513f8f94c67ec6d1c8b6d5bede59b81f6da8c2081f57e49fb28d9784.
Run controls on the workstation before attributing zeros to library chemistry.
Passing controls still does not prove a raw-library census or absence of ionizable
groups in deferred/unsupported molecules. Preserve the original environment.

## Read-only diagnostic

New script: scripts/review_e096_results.py. It reads small reports, blocks.csv,
transform.json and split-decisions.json, verifies their seals and checks population,
size bounds, child fractions, feature summaries, axis participation and E095 size
comparison. It never opens the large SQLite files or MOL2 library. Optional four
synthetic controls exercise the installed factory and residue mapping. Use a fresh
output outside the frozen E096 folder. The diagnostic does not repeat independent
full membership validation or compare individual E095/E096 memberships.

To avoid modifying the running-code checkout, fetch and extract only the diagnostic:

```bash
git fetch origin feature/structure-guided-chat
git show origin/feature/structure-guided-chat:scripts/review_e096_results.py > /tmp/review_e096_results.py
BASE=/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000
PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}" python /tmp/review_e096_results.py \
  --root "$BASE/final-blocks-e096-joint" \
  --baseline "$BASE/final-blocks-e095-properties" \
  --feature-self-test \
  --output "$BASE/e096-review-$(date +%Y%m%d-%H%M%S).json"
```

Run from the Group-Agent checkout in the original aidd-workstation environment.
Return the small diagnostic JSON. If a self-test fails, inspect the original feature
definition and sampled raw structures before modifying frozen descriptors. If it
passes, raw structural inventory is still needed to explain both zero populations.

Next scientific step after diagnostics: fixed molecule-held-out evaluation of E095
and E096 retrieval behavior, then PLANTS/MDM2 and user N-E rescoring on representative
samples. Preserve special-pool search. No sampled result authorizes rejection of
unsampled blocks. Binding-affinity model is still undecided.

Local validation: 16 targeted spatial/diagnostic tests passed; English guard passed.
No frozen E096 production module was modified. MOLIQ and workstation outputs remain
unchanged. This review records user evidence separately from local fixture checks.
