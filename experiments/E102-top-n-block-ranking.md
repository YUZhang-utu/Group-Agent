# E102: molecule-level Top-N block prioritization

Date: 2026-10-05. User explicitly confirmed distinct-molecule Top-10 means.
This is a requested analysis implementation, not a new docking experiment.

For each scheme/block/receptor, group sampled conformers by source molecule ID;
take the lowest ChemPLP score per molecule, select the ten lowest molecule scores,
and use their arithmetic mean to order blocks within that scheme/receptor.
Also report Top-5 and requested-N means. Stable identity ordering resolves ties.
Whole-block means remain historical descriptive outputs, never ranking keys.

The primary ranking requires N successful units and a complete panel. Fewer
units yield a null Top-N mean, not a smaller denominator. Missing scores are
retained and incomplete panels have no primary rank. Molecules crossing blocks
remain dependent. No scheme superiority, experimental affinity or enrichment
claim is made from these exploratory priorities.

Verification design: a synthetic block with a favorable top tail but an adverse
whole mean must outrank a uniformly moderate block; duplicate conformers cannot
occupy multiple molecule slots; the best conformer must be retained even when
its molecule temporarily leaves the bounded candidate set; failed/undersized
panels cannot receive competitive ranks. Test deterministic ties, custom N,
artifact tampering, session ownership, Chat query routing and engine non-dispatch.

Implementation: block_ranking.py; completion-adoption and standard Chat analysis
both append ranking artifacts without changing block_plants.py source or old
docking seals. New bounded ranks/top tools and provider-independent slash commands
expose the saved candidates. Focused engineering regression: 56 passed. Full
regression result is recorded in research_log.md after completion. No scientific
model inference or actual workstation ranking has been run here.
