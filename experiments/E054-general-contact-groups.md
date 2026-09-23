# E054: General contact evidence and same-pose grouped scoring

Protocol written before implementation and validation, 2026-09-23.

Hypothesis: separating exhaustive quality-aware atom proximity evidence from
sparse pharmacophore features prevents unrepresented contacts from being called
absent. Explicit evidence-bound spatial groups and optional-family maxima avoid
residue overlap and duplicate-mode scoring without requiring all observed contacts.

Implementation contract:
- Target identity, canonical mapping, chain, model, atom identity, coordinate
  frame, source hashes and uncertainty accompany contact evidence. No MDM2
  residue names or pocket groups in the general algorithms.
- Preserve all nearby atom pairs and all donor/acceptor proximity alternatives.
  Proximity is not a validated hydrogen/halogen bond, affinity or induced fit.
- Quality-excluded atoms remain visible as unknown. Unsupported chemistry is
  reported rather than guessed. Old sealed evidence stays unchanged; rebuild or
  augment into fresh artifacts.
- Spatial group geometry comes from explicit atom selections in verified source
  complexes in the reference frame, not LLM-invented coordinates. A candidate
  atom has at most one spatial group; ambiguous boundary atoms are unassigned.
- Optional mode families use a maximum, count once and are disjoint from other
  optional families/individual weights. Hard OR groups retain their old meaning.
- Nonlinear occupancy rewards are explicit, monotone, configurable by group
  count and evaluated on one pose. There is no default hard occupancy minimum.
- Legacy designs retain identical scoring when extensions are absent. All new
  fields must pass local validation through chat, adoption and execution.

Verification: synthetic non-MDM2 atom/residue identities; chain/alternate/model
ambiguity; missing atoms; chemical class limitations; common rigid transforms;
overlapping regions; same-pose rather than molecule-union occupancy; duplicate
family rejection; legacy score equivalence; malformed profile rejection; real
48-complex MDM2 contact audit and explicit inspection of the three spatial
exceptions. Check new reference control scores before any library run.

No library scan or live LLM call is required for this development. Cross-target
fixtures verify software generality, not prospective biological validity for all
targets. Record unsupported cases and exploratory results explicitly.

Follow-up protocol after the initial coordinate audit: recheck the 48 existing
instances under a severe reference-state overlap gate (the existing screening
cutoff of 1.2 A), retain incompatible cases as pending, and rebuild consensus into
fresh outputs from cached raw evidence with networking disabled. Compare a design
containing the incompatible 4JV9 template with a reference-only design. The former
must require template/state review; the latter must remain executable. This is
not a universal steric model or a reason to label the pending ligands inactive.
