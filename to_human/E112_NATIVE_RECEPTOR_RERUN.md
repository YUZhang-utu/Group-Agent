# E112 native-receptor rerun

Use the activated aidd-workstation environment, from the repository root. Original workstation coordinates, library artifacts and E109/E111 outputs are required; a reports-only archive is insufficient.

## Fresh bounded pose search (recommended)

    git pull --ff-only origin feature/structure-guided-chat
    python scripts/audit_native_state_poses.py --previous "$HOME/mdm2-e111-strict-2000" --prepared "$HOME/mdm2-e109/prepared" --reference "$HOME/mdm2-e109/pilot-full" --output "$HOME/mdm2-e112-native-fresh-200" --fresh --sample 200 --workers 8 --seconds 600

Adjust input paths only if the original directories differ. Use a new output directory on every run. Python must come from the already activated aidd-workstation environment; do not override it with dl4s.

Each reference ligand and its own target-chain receptor undergo the same verified rigid transform. Native contact hypotheses are rebuilt using the adopted residue/atom/type inventory. Every accepted seed must simultaneously satisfy Ile61 CD1 hydrophobic contact, at least three distinct residues, shape >=0.50, feature color >=0.30, minimum heavy-atom distance >=2.0 A, contact-envelope coverage >=0.65 and each of three geometric groups >=0.50. No thresholds are relaxed.

First 200 entries of the pre-randomized E111 conformer sample are searched using bounded feature-pair/principal-axis seeds. This generates new orientations for existing conformers; it is not docking, torsion generation or energy minimization. Original coarse eligibility is preserved. Eight processes share a cooperative 600-second search budget; source hashing, native setup and in-flight calculations add overhead. This is not a measured total-runtime guarantee.

States lacking the specified native contacts are explicitly unsupported. If all states are unsupported, the run stops with diagnostic controls and an incomplete report. Timeout produces no overall pass fraction; partial conformers are excluded from complete counts. A complete result may legitimately have zero hits. Rates describe sampled conformers, not complete molecules.

Native reference controls and ordered cross-template shape/color comparisons are saved. Self-similarity is not evidence that thresholds recover diverse active ligands. Contact matching remains geometric; sterics use an untyped distance floor. Coverage groups are a reference contact-envelope proxy, not canonical subpocket or cavity-volume measurements. Receptor scope is the verified target chain, without extra symmetry mates or waters.

## Optional paired diagnosis

Omit --fresh and use a different output directory to reevaluate the exact E111 saved poses in native receptors. Their shape/color scores remain unchanged. This answers a different question from fresh search.

## Return results

    tar -czf "$HOME/e112-native-report.tar.gz" -C "$HOME/mdm2-e112-native-fresh-200" .

Inspect report.json, reference_controls.json, native_states.json, cross_template_scores.csv, conformers.csv and passing_poses.json. Receipts preserve source and code hashes. Original scientific modules and input files are unchanged.

## Validation status

32 synthetic and related regression tests passed locally; English-content guard passed. Real native structures and library artifacts are unavailable on this Windows host. Workstation execution and measured pass rate remain pending.
