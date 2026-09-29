# E093 receptor selection advice

The existing agent `pockets` step now delivers a short receptor recommendation,
not a raw cluster table. It reports qualified chains, independent PDB entries,
generic scaffold-series proxies, the proposed receptor set, coverage alternatives,
uncovered structures and explicit limitations. The user retains the selection.

## Decisions

- Limited structures or classified scaffold series: provisional single reference.
- All observed differences within an adequately sampled near-repeat proxy:
  provisional single receptor. This is not evidence of rigidity.
- A nontrivial partition with supported members, a cutoff plateau and exact
  leave-PDB-out stability: partition representatives, not proven physical states.
- Otherwise: PDB-weighted coverage representatives and a recommended count.

Global and smooth local geometry are blended equally after zero-preserving
`x/(x+s)` normalization, with `s=max(0.1,p90)` fitted separately and frozen in the
report. Chemical weight remains 0.15. The old max metric and partitions remain
diagnostics. Neither weights nor tolerances are validated target-independent
physical thresholds. Grid/alignment perturbation calibration remains pending.

Near-repeat evidence uses same-entry chains, with one p90 contribution per entry.
At least three repeat entries are required to use this provisional tolerance;
otherwise an explicitly uncalibrated engineering default is used. Coverage is
not recall, equilibrium population or ligand accessibility validation.
Series are generic Murcko ring-scaffold proxies; missing/acyclic ligands remain
unclassified. Polymer chemistry is not silently inferred.

## Agent and adoption

Use the existing `pockets` intent after structure diversity. The `reference`
object accepts `maximum_representatives` (default 32) and `coverage_fraction`
(default 0.95), in addition to the explicit reference ligand and chain.
The budget is a ceiling, not the claimed optimal number of receptors. A shortfall
is displayed. The normal adoption operation accepts the recommended set; actual
IDs from `receptor_options` allow individual structures to be chosen instead.
No receptor is adopted automatically. Viewer selection follows adopted IDs.
Coverage consensus uses overlapping local neighborhoods, not forced state bins
or a global intersection of all receptor contacts. Representatives without a
prepared ligand reference are marked as requiring preparation before consensus.

The domain agent reads the recommendation, uses its existing Europe PMC
`literature_search` tool for target-specific support, cites returned sources and
separates publications, computed observations and hypotheses. Offline reports
explicitly say literature is not reviewed. Live provider behavior remains to be
accepted on the workstation; prompt wiring is not a completed literature review.

## Existing corrected reports

Re-evaluate without recomputing pocket coordinates or touching an old report:

```bash
python scripts/review_receptor_advice.py \
  --source /path/to/corrected-v2/report.json \
  --output /path/to/new-receptor-advice \
  --maximum-representatives 32 --coverage-fraction 0.95
```

All original source hashes are checked. The output contains a new sealed matrix,
`report.json`, `report.html`, and `receptor-recommendation.txt`. A new normal
`pockets` task uses the same recommendation code automatically. Restart Chat
after updating code. Existing adopted reports remain unchanged.

## Exploratory MDM2 replay

67 chains, 59 PDB entries, 25 scaffold-series proxies. No robust partition was
resolved by the implemented criteria. Current provisional tolerance yields
approximately 72.3% coverage with 8 references, and 95.2% with 21 references.
The latter is a coverage proposal, not 21 biological states. Five structures
remain outside the selected neighborhoods; they are retained for user selection.
Several coverage references require ligand preparation before interaction
consensus. No docking, live literature review or workstation acceptance is
claimed by this replay. The slow E091 library job is untouched.
