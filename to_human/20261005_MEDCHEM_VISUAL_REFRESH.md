# MEDCHEM visual revision - 2026-10-05

This local UI revision supersedes the first E098 masthead. The reference image is
D:/agent/refer.png; the supplied icon D:/agent/med.png remains unchanged.

## Design

- Remove the hero description, structure captions, redundant introductory text
  and persistent right-hand task column.
- Center an italic MEDCHEM wordmark. A real, RDKit-valid cyclic-peptide chemical
  graph, with artistically constrained D-shaped 2D coordinates, replaces D.
- Use a schematic aqueous upper field with water bonds/formulas and a warm nonpolar
  lower field with hydrocarbon skeletons/formulas. It is a design motif, not a density,
  partitioning or solvation simulation. There is no phospholipid illustration.
- Scatter six stable illustrative noncanonical cyclic peptides around the wordmark.
  Examples include N-methyl, fluoroaryl, thienyl, cyclopropyl, pipecolic, proline,
  alkyne, beta-amino-acid, pyridyl, methoxy and biphenyl substitutions. They are valid
  graphs, not identified experimental compounds, library hits or synthesis claims.
  Exact SMILES and provenance descriptions remain in web/peptides.json, without
  visible captions. scripts/build_chat_molecules.py regenerates the static SVGs.
- Tasks & results is collapsed by default. Mouse hover previews the floating panel,
  click pins it, and Escape/close/outside click dismiss it. Touch uses tap. Existing
  attachments, job controls, viewer artifacts and result cards remain in this panel.

## Validation

Credential-free desktop (1440x1100) and mobile (390x844) browser checks cover images,
minimal masthead copy, drawer hover transfer/pinning/closing, keyboard focus/Escape,
touch, existing conversation submission/persistence, attachment type switching,
workflow coverage and no horizontal overflow or JavaScript errors. Screenshots are
in data/e099-ui. Ten existing Chat/API tests and JavaScript syntax checks pass.
The E098 full scientific regression is historical (743 passed, 6 skipped); this
visual-only revision does not rerun scientific experiments or change docking logic.

## Apply to the existing workstation checkout

Transfer and extract MEDCHEM_E099_UI_UPDATE.zip into a separate directory. It is
cumulative: its hash-checked installer accepts the original source or the previously
delivered E098 update, backs up replaced source and refuses other local edits.
Stop the Chat server after checking that no scientific task is active, then run:

```bash
python install_medchem.py --repo /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
python install_medchem.py --repo /mnt/medchem_taltio/wrk/yu_agent/Group-Agent --apply
```

Restart the same Chat storage/session using the command in
20261005_MEDCHEM_CHAT_HANDOFF.md. The actual 8,445-job / 821,619-score docking output
is unchanged. No remote deployment, Git publication or new docking took place here.
The completed-result attachment and prompt-analysis instructions remain the same.

## Typography follow-up

The E100 local revision keeps the entire illustrated wordmark unchanged and uses
Helvetica, Helvetica Neue, Arial, Liberation Sans, sans-serif for the interface.
Helvetica is requested from installed system fonts; no proprietary font is bundled.
Body text is 15px/500, navigation and controls are 13-14px with emphasized controls
at 600, headings are 18-19px/600, and UI text is black or near-black. Technical paths
and status labels also use the same family. Light primary button fills keep black
labels readable. Mobile text fields are 16px. Desktop/mobile interaction and layout
checks pass; screenshots are in data/e100-ui. The scientific backend is unchanged.
The latest cumulative package is MEDCHEM_E100_TYPOGRAPHY_UPDATE.zip; it accepts the
original, E098 or E099 delivered source hashes. Use install_medchem.py as above.
