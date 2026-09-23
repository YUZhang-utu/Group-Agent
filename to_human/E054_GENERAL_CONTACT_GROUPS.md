# General contact evidence, receptor states and grouped scoring

Implemented on feature/structure-guided-chat, 2026-09-23.
The general algorithms contain no MDM2 residues, subpocket names or target IDs.
MDM2-specific definitions live only in an explicit regression script and reports.

## What changed

1. Every new consensus includes a sealed `contacts.json` ledger. It records all
   quality-observed ligand/target heavy-atom pairs within 4.5 A, original atom
   identities, canonical residue mapping, the common frame, CCD-perceived ligand
   chemistry, excluded/uncertain atoms and source hashes. Atom pairs are not
   restricted to existing sparse pharmacophore features or the nearest polar
   partner. Raw proximity, carbon/sulfur proximity, halogen proximity and polar
   donor/acceptor proximity remain hypotheses, not validated interaction energies.
2. Admission now checks severe ligand overlap against the fixed reference
   receptor at 1.2 A, matching the existing screening overlap cutoff. Conflicts
   stay pending with evidence and reasons instead of joining the same-state
   cohort. This is not complete VDW validation, alternate-state preparation or
   an inactivity label.
3. Optional mode families use their best member once. Their weights share the
   optional-score normalization with individual optional anchors. A member cannot
   also appear in another family, an individual optional weight or a hard rule.
4. Spatial regions use explicit ligand atom selections from compatible, sealed
   source complexes. The server resolves coordinates; arbitrary points are not
   accepted as design edits. Each candidate atom is assigned to the nearest
   qualifying region, with near-ties left ambiguous, and counts at most once.
5. Explicit nonlinear occupancy rewards are added to the same-pose score. There
   is no implicit all-regions hard requirement. Existing anchor, pocket and pose
   cutoff rules still apply: these scoring extensions do not remove the existing
   requirement for at least one matching selected anchor.
6. Both per-seed and molecule/template aggregation preserve distinct spatial
   occupancy signatures even when their anchor masks are identical. Payloads
   remain individual poses; a molecule-level union is never a simultaneous pose.
7. Selected templates must pass crystal self-controls. Otherwise adoption returns
   `needs_template_state_review`; the funnel refuses that design. A completed
   design report alone no longer implies readiness for screening.

## Editable design fields

Old designs remain valid with unchanged numeric scoring when extensions are
absent. New fields are optional but locally validated in all routes:

| Field | Meaning |
|---|---|
| optional_groups | List of `{id, anchor_ids, weight}`; best member counts once |
| spatial_groups | List of `{id, reference_query, ligand_atoms, radius, minimum_atoms}` |
| occupancy_rewards | One monotone value in [0,1] per occupied-group count, beginning with zero |
| occupancy_weight | Relative contribution of the occupancy reward, in [0,1] |
| spatial_ambiguity | Distance difference in A below which competing qualifying regions are ambiguous; default 0.5 |

All reference queries and atom names must be present and quality-observed in the
sealed ledger and compatible with its target/reference state. A reference atom
cannot define two groups. Region radii are explicit, between 0.1 and 4 A; minimum
atom counts are explicit positive integers. The current limits are 16 spatial
groups and 20 selected feature anchors including optional family members.

For three regions, `[0, 0.1, 0.3, 1.0]` is an example nonlinear reward, not a
validated universal default. Rewarding three regions more strongly than two
does not guarantee strict lexicographic ordering over all other score terms.
Weights normalize as:

`score = (gaussian_weight * gaussian + optional_weight * optional + occupancy_weight * occupancy_reward) / sum(weights) - exclusion_penalty`

Adding occupancy weight 0.3 to the old 0.7/0.3 weights changes their normalized
contributions; it does not preserve Gaussian at 70 percent of the new total.
No production coefficients or MDM2 design have been automatically adopted.

Spatial occupancy here is an explicit geometric overlap measure, not proof of a
specific hydrophobic interaction. Atom chemistry and polar direction remain
separate evidence. The current region references are admitted nonpolymer ligand
atom selections, not automatically discovered cavity boundaries or peptide sites.

## Validation and the three reviewed spatial exceptions

- Full suite: 412 passed, two skipped. Includes non-MDM2 synthetic identities,
  arbitrary chains/residue numbers, multiple polar partners, halogen evidence,
  source/frame tampering, uncertain coordinates, grouped scoring, rigid-transform
  invariance, boundary ambiguity, hard-vs-optional rules, persistence and legacy
  score equivalence. These establish software behavior, not biological validity
  for every target.
- Historical cohort: 48 complexes / 45 PDBs, 6905 ledger pairs exactly match the
  independent atom audit and distances. Includes 3671 carbon/sulfur, 1153 halogen
  and 59 donor/acceptor proximity labels; categories are not mutually exclusive.
- 4JV7 and 4JV9 are about 2.70 A from the nearest p53 Leu26 probe atom. More
  decisively, their aligned ligands collide with fixed-reference GLN59 atoms at
  about 0.64-0.75 A. They are pending receptor-state cases, not inactive molecules.
- 4ZFI's closest Phe19-probe candidate atom is C7 (2.249 A), but it is closer to
  the Trp23 probe (1.623 A). It covers three regions when using reviewed 5C5A
  ligand atom selections instead of p53 side chains. Its apparent exception is
  reference-dependent; it remains admitted and must not be rejected for this.
- 4MDN was additionally found incompatible with the fixed TYR100 state (minimum
  overlap about 0.474 A). It is also pending.
- Offline consensus rebuild: 45 prepared complexes, 119 feature modes, 39 unique
  templates, 6507 atom pairs, zero preparation failures. A subsequent ledger
  augmentation confirms zero unresolved reference-state checks in this cohort.
- A deterministic fixture with grouped bonuses and two templates correctly
  blocks incompatible 4JV9. Its reference-only counterpart is ready. A fresh
  45-complex-cohort fixture with 20 frequency-selected optional anchors, one
  optional family and reviewed spatial groups passes 45/45 crystal controls.
  Its reference-derived cutoff is about 0.290259, an exploratory fixture value.

These fixtures are neither live LLM recommendations nor user-adopted designs.
No full-library scan, full-library NH count, new independent-active recall,
docking, affinity validation or induced-fit prediction was performed.

Local receipts (ignored raw artifacts, retained in the development workspace):

- `data/e054-contact-ledger-v3/`: historical cohort with state diagnostics.
- `data/e054-consensus-rebuild/`: full offline rebuild from raw cached structures.
- `data/e054-compatible-contact-ledger/`: current compatible cohort and full ledger.
- `data/e054-groups-validation-final/`: independent pair comparison, two-template
  failure control, compatible/fresh adoption fixtures and per-pose spatial counts.

Reproducible entry points are `scripts/rebuild_cached_consensus.py`,
`python -m aidd_agent.contact_evidence`, and
`scripts/validate_mdm2_contact_groups.py`. The latter is explicitly an MDM2
regression case; it is not the general engine or a new scientific recommendation.

## Workstation continuation

Update the existing feature branch and restart Chat with its existing storage,
runtime and provider profiles. Use a fresh task after the code update; old pose
databases are not migrated or resumed under the new schema.

Reuse completed structure diversity. In the same conversation, rebuild consensus
with explicit target/reference identity, for example:

> Rebuild the pocket consensus for human MDM2 Q00987 using reference ligand
> 5C5A:NUT:A:201 and target protein author chain A. Include the new contact ledger
> and reference-state checks. Do not start library screening.

Inspect the ledger, pending states and new templates, then request a fresh
`/recommend`. Rebuilding changes the cohort and may change anchor IDs: do not
blindly reuse the previous 20 IDs or a proposal sealed against older evidence.
Ask for edits using current IDs and explicit reviewed atom selections/rewards.
Edited designs are created directly; `/adopt` without edits accepts the original
proposal, so do not issue it afterward expecting to preserve a separate edit.

The conditional LEU54 mandatory proposal remains deferred pending actual library
NH/feature retention and independent-positive checks. GLN72 O and OE1 remain
distinct. VAL93 polar proximity is not automatically promoted to a validated
hydrogen-bond rule. Unsupported chemistry, missing/alternate atoms, water/metal
geometry and alternative receptor states retain explicit limitations.
