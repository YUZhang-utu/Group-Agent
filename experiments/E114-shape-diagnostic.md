# E114 bounded shape diagnostic

2026-10-20 exploratory protocol, before analysis. Question: does the zero-hit E113 sample reflect whole-ligand shape size penalties or insufficient local rigid refinement, rather than established absence of suitable molecules?

Use verified E113 receipt only for local analysis. Recover nine reference heavy-atom point clouds from self-control recovered poses with shape/color close to one. Validate reconstructed reference shapes against saved shape scores for every analyzed pose; abort on mismatch. No atom typing or receptor reconstruction is inferred.

Deterministic panel: all 35 geometry-qualified near-miss conformers (their best coverage-with-geometry pose), six leave-self-out control examples chosen by best shape where present, plus the top 20 library best-shape poses by score (deduplicated by task/pose). Include all available leave-self-out cases rather than select outcomes. This is an outcome-selected diagnostic panel, not a new random prevalence estimate.

Report atom count ratio, original Gaussian self/cross overlaps, self-norm Cauchy-Schwarz upper bound on Tanimoto, asymmetric point-neighborhood coverage both directions at 1.5 A, and reference contact-region coverage. Cauchy bound uses current point-cloud kernel, not a volume or affinity bound. Refined coordinate kernels use sigma 1 A and no cutoff, unchanged from E113.

Perform one bounded local rigid refinement from each saved aligned pose, with rotation-vector components within +/-0.5 radians and translation components within +/-2 A about candidate centroid, max 180 objective evaluations per case, total cooperative 120-second budget. No torsion generation, global alignment or receptor energy optimization. Report optimizer success/status and objective calls; retain the better original/refined shape. Deadline checked during objectives. Partial rows/status explicit.

Local optimized poses are geometry-unvalidated: receptor clash and contacts cannot be rechecked without workstation artifacts. Export composed original-candidate-to-reference transforms for library cases only. A short optional workstation checker should re-evaluate these proposals against sealed native receptor and chemical artifacts, checking contacts, whole color, geometry and shape. No library-wide seed search. Keep original thresholds unchanged and do not report optimized shape crossings as hits before that check.

Also compare kernel width 0.75/1.0/1.5 A on unchanged poses diagnostically, without transplanting the 0.50 cutoff across widths. No threshold selection based on this panel. No conclusion that there are no suitable molecules in all 10000 or the full library: the run only covers finite conformers/seeds and adopted templates.
