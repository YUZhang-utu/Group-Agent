# E082 primary source notes

Accessed 2026-09-28. These support design choices, not project validation results.

- Bista et al. (2013), "Transient Protein States in Designing Inhibitors of the
  MDM2-p53 Interaction", Structure 21:2143-2151.
  https://doi.org/10.1016/j.str.2013.09.006
  Primary structural work reports ligand-associated Tyr100 and lid-region states.
  Project implication: inspect receptor pocket states instead of relying on one
  ligand-bound geometry. Do not generalize any one state to all MDM2 ligands.
- "Small Molecule Inhibitors of the MDM2-p53 Interaction Discovered by
  Ensemble-Based Receptor Models" (2007), JACS.
  https://doi.org/10.1021/ja073687x
  Reports multiple receptor models and experimentally identified distinct inhibitor
  scaffolds. This supports testing the ensemble hypothesis, not assuming our
  library or PLANTS/N-E combination has the same performance.
- Korb, Stutzle and Exner (2009), "Empirical scoring functions for advanced
  protein-ligand docking with PLANTS", JCIM 49:84-96.
  https://doi.org/10.1021/ci800298z
  Primary description of PLANTS empirical scoring. A docking score should retain
  its engine-specific semantics rather than be relabeled measured affinity.
- AutoDock Vina official macrocycle documentation:
  https://autodock-vina.readthedocs.io/en/stable/docking_macrocycle.html
  Demonstrates that some docking procedures change ring geometry. This motivates
  source/output pose attribution; it does not establish PLANTS ring capabilities.

Recommendations about fixed receptor budgets, unfiltered sampling, same-pose
scorecards and held-out cross-docking are project protocol inferences. Validate
them empirically before deploying live block rejection.
