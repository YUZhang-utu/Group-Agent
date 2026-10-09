# E107: independent block-search benchmark

Protocol 2026-10-09. Retrospective/exploratory: the twelve arm summary counts were
already inspected. This is not a preregistered confirmatory experiment.

Question: which partition and block-priority scorer provides useful coverage at
a fixed retrieval budget? Preserve E106 thresholds, query templates, identities,
Top-10-distinct-molecule block statistics and complete-block stopping. No new
search, docking or rescoring is launched by this protocol's initial audit.

Primary comparison: six Top-10-block-limit arms, E094/E095/E096 crossed with
ChemPLP/EquiScore. Each exports 100000 conformers. Top-5 is the routing-limit
ablation, not a change to the ten-molecule block statistic. Shortfall arms remain
shortfalls; no padding or outside-block candidates. A secondary common-size
prefix audit can use 76933 conformers across all twelve arms, after verifying
candidate file order against the export ranking implementation.

Stage 1: verify complete reports, selected candidate file hashes, unique CIDs,
CID-to-molecule consistency and reported counts. Report export completion,
survival fraction, searched workload, molecule counts and conformers per molecule.
Compute conformer/molecule intersections, Jaccard and directional containment;
test exact Top-5/Top-10 set equality. No biological superiority follows.

Stage 2: join verified original molecular structures and partition metadata.
Report building-block, stereochemistry and cis/trans state coverage separately.
Never collapse molecules by formula. Morgan radius 2, 2048 bits, chirality enabled
is a proposed 2D descriptor, not a measure of conformational cis/trans coverage.
Compare distributions at matched molecule sample sizes with a fixed recorded seed;
avoid a quadratic all-pairs matrix across hundreds of thousands of candidates.
Missing/unsupported chemistry must be reported, not silently removed.

Stage 3: reference retrieval benchmark requires a whole-library run with identical
query, thresholds, engine, ranking and library eligibility. Freeze reference Top-K
before inspecting arms. Report reference-conformer and reference-molecule recall,
explicitly as retrieval recall, not active-compound recall. Audit whether RRF ranks
are block-local; do not assume block-local ranks equal global-library ordering.
Include no-block full-library and repeated random-block/uniform-candidate controls
at matched evaluated work or cost. Uniform molecule and conformer baselines are
different. Report the regular-block universe separately from excluded special pools.

Stage 4: experimental labels must be independent of query design and block scoring.
Confirm labeled molecules exist in the evaluated library, preserve stereochemistry,
define assay-compatible activity thresholds, and separate unknown from inactive.
On a fixed labeled universe, compute molecule-level active recall, precision and
EF at fixed selected budgets. For EF1%, use the same denominator N and
ceil(0.01*N) across methods, not each method's selected subset. Rank-based BEDROC
requires a declared whole-universe ranking/missing-rank policy. Query templates
and molecules used to determine priorities are leakage controls, not held-out hits.
Bootstrap paired molecules, retaining all their conformers together; repeated
targets/routing samples are needed for robustness beyond the current single run.

Timing: report cache-hit and cold timings separately. All current Top-10 arms have
cache reuse; do not rank their wall times as fresh algorithm performance. The two
shortfall Top-5 arms also cannot be compared as completed 100000-candidate runs.
For deployment report amortized cost; for methods report partition/index build,
initial block sampling, PLANTS/EquiScore prioritization and retrieval cost separately.
Do not silently erase the larger block-sampling cost of finer partitions.

Decision: no arbitrary combined score. Report a Pareto comparison of independent
reference/active recovery, cost and coverage; select a policy only after choosing
the intended application budget and validating on held-out targets.

Implementation sources:
- https://rdkit.org/docs/source/rdkit.ML.Scoring.Scoring.html
- https://www.rdkit.org/docs/source/rdkit.Chem.rdFingerprintGenerator.html

Current prerequisite: remote candidates and baseline/label availability pending.
User-provided receipts confirm all twelve E106 jobs complete; local Windows has
no direct workstation connection. A portable read-only audit is supplied separately.
