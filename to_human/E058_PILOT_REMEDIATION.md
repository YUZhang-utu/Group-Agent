# E058 pilot diagnosis and equivalent acceleration

Run from the repository in the workstation environment:

```bash
bash scripts/run_e058_review.sh
```

The script reads the existing 256-molecule pilot, recovers up to five failing
original MOL2 records, replays saved poses against both region definitions, then
reruns the original 256-molecule seed/sample settings with compiled assignment.
It never runs the whole library or changes screening thresholds. Requires Numba
for the final rerun; if absent the two preceding diagnostic reports remain valid.
Install it in the active environment using `python -m pip install numba` only if
needed. The compiled backend is explicit; missing compilation support is an error,
not silent fallback. The default production backend remains Python.

## Chemistry

The original MOL2 loader permits a separately recorded aromatic-no-kekulize mode.
Scaffold reconstruction applies strict sanitization to a heavy-only graph without
stored hydrogen counts. These are two distinct candidate causes. Source recovery
verifies registry record content hashes, heavy atom order/numbers, charges,
aromatic flags and bond orders. It reports aromatic heteroatom hydrogen counts.
Neither mode proves incorrect input chemistry nor licenses guessing missing H.
Actual failing workstation records remain pending until this command runs.
Scaffold failures are preserved and now counted separately from valid groups.

## Selectivity

The baseline has 226/256 passing molecules. Template exclusive hits quantify
which templates add unique members to the union, independently of template order.
Score/contact-count summaries cover saved representatives only. They cannot
identify all rejected poses or predict relaxed-threshold results. Current adaptive
coarse bounds use the receptor-state cohort, and final membership requires at least
one selected contact plus the composite threshold. No occupancy rule was active.
The Gaussian stage label now explicitly says evaluated, not passed.

## Two definitions

The portable MDM2 review fixture preserves prior E053 peptide alignment and E054
5C5A atom selections. It uses a 2 A radius, one assigned atom and 0.5 A ambiguity
margin as explicit exploratory settings. Atoms compete only between eligible
regions, matching current contact_groups semantics. This detail matters near
boundaries and is not interchangeable with competition against every region.
Both definitions occupied means agreement, only one means boundary, neither means
unoccupied, missing means unknown. Agreement is not independent confidence.
These are geometric occupancy definitions, not proof of hydrophobic contacts.

Local aligned crystal controls: 44/45 have three agreeing regions. 4ZFI has two
agreeing regions and boundary C_Phe19 (p53 false, 5C5A true). No hard three-region
gate or numerical boundary bonus is adopted. Workstation replay only sees saved
passing representatives; occupancy-aware ranking could select additional poses.
Thus replay is a review tool, not an exact alternative production selection.

## Performance and acceptance

An optional Numba compilation of the existing Hungarian loop preserves operation
order and ties without fastmath. Local randomized/tied/empty rectangular cases
match exactly. A 300-case 14-by-30 matrix microbenchmark measured 0.162 s Python
versus 0.00179 s compiled, excluding compilation: not an end-to-end speed claim.
Seed generation is unchanged and remains a significant bottleneck.

Compare reports require identical sample, design and complete pose payloads.
Compilation warmup is recorded. Only an equivalent complete workstation rerun
supports acceleration claims; incomplete samples cannot be extrapolated.
After source recovery and saved-pose occupancy review, decide chemistry repair
and region scoring before any larger run. Biological recall remains deferred.

Validation: 452 tests passed, 2 skipped in the final stable-source full run;
25 explicit compiled-backend preselection/consensus/remediation tests passed.
