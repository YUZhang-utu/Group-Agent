# Ideas and open questions

## Active E019 questions

- For 60-dimensional USRCAT, which of m20, m30, or zero-padded d64/m32 gives the
  best recall/build-time/resident-memory tradeoff?
- What nprobe and L1 top-N reach acceptable exhaustive top-1000 recall at one
  million and at larger physical shards?
- Can FAISS IDs be represented implicitly/customly without losing resumability,
  or must the local RAM budget reserve roughly 8 GB per billion conformers?
- Which partial RDKit sanitization operations are necessary and sufficient for
  feature equivalence to the fully sanitized reference path?
- Are source conformers contiguous under a trustworthy source molecule ID, and
  what disk-backed lookup is needed if they are not?
- Does post-z-score HBD/HBA block weighting improve candidate recall without
  reducing scaffold diversity?
- Which anchor definition is most reproducible across protonation, water-
  mediated interactions, and uncertain crystal density?
- Do projected hydrogen-bond points improve EF1% over atom-centered color points
  after controlling for the same candidate budget?
- Does query-biased ColorTversky improve scaffold diversity without admitting
  excessive unmatched polar functionality?
- Can one/two terminal-bond relaxation recover sparse-ensemble misses without
  turning similarity refinement into docking?
- What approximation order provides the best ranking/runtime tradeoff for
  Gaussian molecular self- and cross-overlap?
- Does the locked `W/n` normalized-mixture representation need an additional
  group-level self-overlap normalization, or is empirical Tversky calibration
  sufficient without improving enrichment?
