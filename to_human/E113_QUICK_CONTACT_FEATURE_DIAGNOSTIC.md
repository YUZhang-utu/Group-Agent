# E113 quick contact-feature diagnostic

Use the activated aidd-workstation environment from the repository root:

    git pull --ff-only origin feature/structure-guided-chat
    python scripts/diagnose_contact_features.py --previous "$HOME/mdm2-e111-strict-10000" --prepared "$HOME/mdm2-e109/prepared" --reference "$HOME/mdm2-e109/pilot-full" --output "$HOME/mdm2-e113-contact-200" --sample 200 --workers 20 --seconds 120

Use a fresh output directory. No new E111 run is needed. This selects the first 200 entries of the sealed, pre-randomized E111 10000 sample before inspecting outcomes. Native structures and original library artifacts must remain accessible. The Windows reports-only workspace cannot execute the scientific pilot.

## Scope and budget

One shared 120-second cooperative search budget covers 22 control tasks followed by 200 library conformers. Controls are self-template and leave-self-template-out bounded searches for the 11 existing native ligands after a fixed rigid rotation/translation. Only supported target states are searched; unsupported self controls are reported explicitly. Controls bypass library coarse eligibility. Source verification/setup and in-flight work add overhead, so this is not a hard wall-clock limit or a promise of completion. Partial tasks have no overall pass fraction and are excluded from complete counts. No automatic enlargement occurs.

## What changes

The new diagnostic averages directional one-to-one match scores over unique reference feature indices supported by the adopted native contact inventory. Extra candidate features do not enter the denominator. This is selected contact-feature coverage, not a claim of essential binding features. No solvent exposure, buried-polar penalty or direct energetic validation is implemented. Shared residue labels for one reference feature do not increase its coverage weight.

Original shape/color and all original gates are retained as the baseline. Sensitivity levels 0.5/0.7/0.9 compare coverage as a replacement for the color gate while retaining shape >=0.50 and all geometry gates on the SAME seed. These levels are exploratory, not adopted production cutoffs. Contact-plus-all-geometry is also counted without shape/color for diagnosis only. Shape remains a potential complete bottleneck.

pose_diagnostics.json saves best shape, best color, best coverage and best coverage among geometry-qualified poses per task, with simultaneous scores, assignments, native-state identity, transforms and aligned heavy-atom coordinates for completed searches. Partial-task records are explicitly incomplete and may lack coordinates. Never combine different records to infer same-pose acceptance. Self-control success only checks search recovery; leave-self-out references can be chemically related and do not establish independent activity enrichment.

## Return

    tar -czf "$HOME/e113-contact-report.tar.gz" -C "$HOME/mdm2-e113-contact-200" .

Return the archive, including report.json, protocol.json, pose_diagnostics.json, sample.csv and native reference diagnostics. Output hashes are in report.json.

## Local validation

37 diagnostic and related regression tests passed. Tests cover extra-feature denominator invariance versus whole-ligand color, missing-reference penalties, one-to-one assignment, shared-index deduplication, directional matching, transformed-control search recovery, baseline score agreement and timeout accounting. Actual workstation result remains pending. Protocol commit: c160d09.
