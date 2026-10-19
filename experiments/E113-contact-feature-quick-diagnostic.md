# E113 fast contact-feature diagnostic

Exploratory protocol before implementation and workstation execution. Hypothesis: whole-ligand color may penalize candidates that reproduce the adopted native contact features; reference-normalized feature coverage removes unmatched-candidate denominator penalties. This is not a validated activity score.

Reuse E112 native-state verification and unchanged source modules. Search first 200 pre-randomized E111 sample entries, 20 workers, one shared 120-second cooperative search budget including controls. Do not enlarge automatically. Setup/hashing and in-flight seeds add overhead. Partial tasks have no overall rate.

Compute mean directional one-to-one match score over UNIQUE native feature indices supported by the adopted contact inventory. Unit reference weights, unmatched references zero. Additional candidate features do not enter the denominator. This only covers selected native contact modes, not all biological contacts or experimentally essential features. No new solvent exposure or buried-polar assignment is claimed.

On the SAME contact-qualified seed record original shape/color, feature coverage, collision and envelope/region coverage. Original full gates remain the baseline. Diagnostic coverage thresholds 0.5, 0.7, 0.9 are a predeclared sensitivity panel, not accepted production cutoffs. Replacement counts retain original shape and all geometry gates but omit whole-ligand color. Also report contact-plus-all-geometry counts without shape/color for diagnosis only.

Save best shape, color and coverage poses per conformer, plus highest-coverage pose among those satisfying all geometry gates, including coordinates and complete same-pose metrics. No interpretation from mixing independent maxima.

Control tasks use the 11 existing native ligands after a fixed proper rotation/translation, applying identical seed generation. One self-template recovery and one leave-self-template-out search per ligand; only supported target states are searched. Controls bypass library coarse eligibility and are labeled. These are internal diagnostic controls, not independent active-vs-inactive validation; unsupported self states are explicit. Controls first, then library tasks, sharing the budget.

Local tests must establish added remote candidate features do not depress coverage, missing reference features do, one candidate feature cannot satisfy two independent reference indices, shared contact labels do not double-weight a feature, and original baseline scores/gates remain unchanged. Real workstation data unavailable locally; deliver runnable script and record pending execution.
