# E117: prospective KRAS reference-gate diagnostic

Question: do reference feature assignments reject poses which already satisfy
all five native contact regions and allele-specific residue 12 conditions?

Use the identical 1,000 spread conformers as E116. Keep direct contacts, shape,
whole color, sterics and coverage unchanged. Compare strict acceptance with a
diagnostic bypass of reference feature matching only. Count per-anchor assigned
and score-qualified poses, required residue/class failures and distinct-feature
failures. Report counts per conformer and unique molecule, with diagnostic-only
candidate examples. Do not promote bypass hits to production candidates.

This is exploratory. Native recovery and unchanged production behavior are
local checks; full-library policy selection requires workstation results.
