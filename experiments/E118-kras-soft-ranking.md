# E118: exploratory KRAS direct-pool reranking and docking handoff

Authorized endpoint: 1,000,000 distinct library molecule IDs per target, ten
100,000-record MOL2 files per target when enough eligible molecules exist.
Do not fill a shortfall with failed molecules. No docking is performed here.

The returned counter retains target/molecule IDs only, not passing conformer IDs,
transforms or scores. Recover finite poses for conformers belonging to that pool;
do not call this cached-pose rescoring. Original counts and strict runs stay intact.

Retain the original same-pose direct-contact chemistry and finite seed generator.
Apply the original 2.0 A receptor-heavy-atom distance floor independently of shape.
This is a gross-clash proxy, not atom-specific VDW or energy validation.
Reference directional matching, shape, whole-ligand color, reference contact
envelope and its three regions are soft metrics. Coverage is not cavity volume.
Average maximal directional matches over reference residue groups to reduce
reference-point multiplicity bias. This is a diagnostic soft score, not a new
validated hydrogen-bond metric. Candidate feature competition remains explicit.

Use an exploratory equal-block mean of shape, contact-feature score and mean
regional coverage for same-pose ordering; color is reported but excluded from
selection to avoid the previously identified whole-ligand penalty. All components
are bounded [0,1]. This transparent choice is not fitted to activity or claimed
optimal. Retain three best distinct conformers per molecule, with deterministic
ties; select the best for the one-record-per-molecule docking handoff.

Select 80% by score and 20% by round-robin coarse graph-fingerprint buckets among
remaining eligible molecules. Buckets use deterministic locality-sensitive graph
fingerprints; they are diversity proxies, not exact scaffold clusters or chemical
identity deduplication. A larger quota can be applied later without pose rescoring.

Before uncapped work: transformed native positive controls must pass; a spread
10,000-conformer pilot records geometry retention, score distributions, throughput
and explicit same-policy lineage. A pilot cannot establish million-molecule yield
or a guaranteed completion time. Full execution requires that pilot receipt.
Atomic chunk transactions and input/code signatures support safe resume. An
optional cooperative hour limit saves progress and does not silently finalize.

Export the selected conformer from the registry's original source MOL2 with
record SHA256 validation. Preserve coordinates, atom/bond data and charges;
replace only the title with a traceable stable alias. Include pose transforms,
scores, source locators and target membership in manifests. Strict RDKit parse
failures are quarantined, never silently exported as accepted. No protonation,
tautomer enumeration, minimization, affinity inference or stereochemical dedup.

Tests: native controls, collision rejection, deterministic ties/diversity, chunk
atomicity/resume, source hash failure, partition size and preserved source content.
Workstation speed and eligible yield remain unmeasured until returned evidence.
