# E062: seed and pocket overhead

2026-09-23 protocol. User receipt: 64 conformers, 33088 seeds, 29740
clash-rejected, 3348 surviving; seed generation 1.07947s, pocket 0.96277s,
Gaussian 0.20678s, assignment 0.42046s, chunk wall 2.78895s. Configuration:
16 workers, 64 conformers/chunk, numba assignment. Single chunk is not a
representative template-wide sample. Optimize repeated axial rotations and
batch independent KD-tree queries. Preserve all seed IDs/order/transforms and
physical predicates. Compare to pre-change code and scalar pocket reference;
benchmark reproducible synthetic panels, then full regression. No claimed
workstation speedup until measured on the same real panel.

Exploratory results: exact ordered seed payloads matched c0a84bb on ten random
panels (default_rng62, 20 candidate/12 query features, coordinates normal*3,
types integers 1..4), caps 1/64/512. Five warm cap512 repetitions: 1.33x.
Pocket scalar equivalence on 60 atoms/517 translations/2000 receptor points:
1.54x synthetic kernel speedup. Full regression: 474 passed, 2 skipped.
Artifact: data/e062/synthetic.json. No measured workstation end-to-end speedup.
