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

## Future docking, prediction, and laboratory-loop questions

- Which redocking and cross-docking thresholds are sufficiently predictive to
  promote a WEE1 docking protocol beyond pose generation?
- How should receptor-ensemble and pose evidence be aggregated without hiding
  a bad clash, high ligand strain, or an unsupported interaction?
- Which prediction models have prospective calibration and applicability-domain
  coverage for the macrocycle chemical space?
- Does Pareto selection improve potency/property/diversity tradeoffs over a
  fixed weighted score in prospective rounds?
- What is the smallest manually validated plate assay that can serve as the
  first complete sample-to-result automation vertical slice?
- Which vendor-neutral protocol representation can compile deterministically to
  the laboratory's actual liquid handler, robot arm, and reader?
- What QC and induced-failure tests are required before unattended operation?
- Which active-learning objective improves valid information per experiment
  over random, similarity, and expert-selection baselines?
