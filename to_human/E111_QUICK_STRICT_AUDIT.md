# E111 quick strict-pose audit

New read-only diagnostic; original E109/E110 scientific modules remain unchanged. Protocol: experiments/E111-strict-cached-pose-audit.md.

## Rules (all on one saved pose)

- Ile61 anchor pocket:74a6061918ece571 >=0.5, plus at least three distinct contacted residues in total.
- Independent shape Tanimoto >=0.50 and feature-color Tanimoto >=0.30.
- Minimum ligand/receptor heavy-atom distance 2.0 A.
- Cover >=65% of matched-template atoms within 4.5 A of the receptor, at candidate distance <=1.5 A.
- Cover >=50% of each of three deterministic spatial groups of those template contact atoms. Unsupported group geometry fails closed and is reported.

These numerical floors are exploratory hypotheses, not activity-calibrated thresholds. Contact-envelope groups are geometric proxies, NOT validated canonical subpockets or complete cavity volume. Sterics use untyped distance, NOT a full VDW energy model. Inspect reference losses before production adoption.

## Run on workstation

Use the activated aidd-workstation environment and original E109 dependencies. Do not hard-code dl4s or change an environment during an active run. Run from the repository:

    python scripts/audit_strict_cached_poses.py --prepared "$HOME/mdm2-e109/prepared" --reference "$HOME/mdm2-e109/pilot-full" --output "$HOME/mdm2-e111-strict-2000" --sample 2000 --seconds 600

Use a fresh output directory. Original query, schedule.sqlite, scores.sqlite, receipts and library/query coordinate resources are required. A reports-only archive is insufficient.

2000 conformers are uniformly sampled BEFORE looking at outcomes; old failures remain in the denominator. Audit processing is limited cooperatively to 600 seconds. Source hashing/setup and the current conformer can add time; this is not a measured runtime guarantee. Partial output has a null overall pass fraction.

Only saved best-per-template poses are checked. Other valid poses may have been discarded, so this measures cached-pose retention, not definitive stricter-search recall. Molecule counts over sampled conformers are not full molecule-level pass rates.

## Return

    tar -czf "$HOME/e111-strict-reports.tar.gz" -C "$HOME/mdm2-e111-strict-2000" .

report.json contains old/cumulative/standalone counts, complete-sample pass fraction, reference losses and timing. poses.csv contains residue sets and geometry scores; regions.json stores the exact region coordinates. No alignment, docking or full-library run occurs.

## Local evidence

Protocol committed before execution: 6a7a8b0. Saved crystal contact scores retain 41/45 poses and 10/11 templates; 5HMK:62Q:A:1001 fails the new contact rule. Geometry gates were not tested on real crystal coordinates locally.

26 synthetic/related tests passed, including same-pose conjunction, distinct residue counting, independent score floors, geometry, read-only audit I/O, timeout accounting and output hashes. These are not real-library results. Actual workstation pass-rate testing remains pending because raw databases and coordinate resources are not present locally.