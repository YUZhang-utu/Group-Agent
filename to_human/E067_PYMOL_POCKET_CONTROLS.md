# E067: chain editing, mixed colors and pocket contact overview

The user confirmed that the desktop bridge connected on 2026-09-25. Later feedback
identified missing chain removal, mixed colors and a single pocket-review command.
These extensions are locally tested; real workstation acceptance remains pending.

After pulling feature/structure-guided-chat, save any desired PSE, close the
bridge-launched PyMOL window, stop and restart Chat with its original storage and
runtime, then reopen the completed consensus task. Both processes cache Python
code; updating files alone does not update a running viewer.

## English prompts and direct fallback

| Prompt | Direct command |
| --- | --- |
| List protein chains in v001 and indicate which are within 5 angstroms of the ligand. | `/view {"operation":"chains","objects":["v001"]}` |
| Remove protein chain B from v001 in the displayed copy. | `/view {"operation":"remove_chain","objects":["v001"],"chain":"B"}` |
| Color the ligand in v001 by element. | `/view {"operation":"color","objects":["v001"],"target":"ligand","scheme":"element"}` |
| Color protein chains in v001 differently. | `/view {"operation":"color","objects":["v001"],"target":"protein","scheme":"chain"}` |
| Color the protein in v001 with a rainbow gradient. | `/view {"operation":"color","objects":["v001"],"target":"protein","scheme":"rainbow"}` |
| Show and label the complete residues within 5 angstroms of the ligand in v001, using element colors. | `/view {"operation":"pocket_view","objects":["v001"],"radius":5}` |
| Show all available candidate contacts and the 5 angstrom pocket residues for v001. | `/view {"operation":"interaction_overview","objects":["v001"],"radius":5}` |

Use the actual catalogue IDs. Chain B is an example, not a recommended deletion.
Ambiguous "redundant" chain requests list chains first; proximity alone does not
establish biological redundancy. Explicit deletion removes protein atoms only,
preserves ligand/source files, and clears obsolete contact lines. Reopening the
task restores the original displayed structures.

Newly opened ligands use element colors by default. Pocket view uses cyan carbon
for residues and green carbon for ligands; N is blue, O red and S yellow. The
protein cartoon remains visible, while residue sticks/CA labels identify the
requested neighborhood. Whole residues are selected if an atom lies within the
radius. Radius may be 1-12 A; this is a display setting, not a screening criterion.

The overview exports a residue list and separate contact tables. Gray dashes show
heavy-atom proximity within 4.5 A; yellow dashes show PyMOL polar candidates within
3.6 A. Counts overlap and must not be summed. Distance labels are hidden to reduce
clutter; exact distances remain in the tables. This is not comprehensive validated
interaction profiling. Salt bridges, pi stacking, cation-pi, halogen bonds,
water bridges and metal coordination are explicitly reported as not evaluated.

Acceptance: check /pymol_status after each command, inspect element colors and the
5 A residue list, ensure chain deletion affects only the requested protein chain,
and reopen to verify restoration. Compare English routing with direct commands.
The underlying display APIs follow the
[official PyMOL implementation](https://github.com/schrodinger/pymol-open-source/blob/master/modules/pymol/viewing.py).
