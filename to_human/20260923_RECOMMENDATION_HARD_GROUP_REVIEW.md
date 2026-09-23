# Review of the user's new recommendation

The user supplied a completed recommendation with one mandatory ILE61 mode and
four hard alternative groups. Each group is an OR internally, but the groups
are combined by AND along with the mandatory mode. This differs from the user's
requested optional enhancements and nonlinear spatial occupancy rewards.

Local replay against the matching rebuilt consensus (45 prepared crystal identity
poses, score threshold 0.5, selected 20 anchors) gives:

| Requirement | Passing crystal poses |
|---|---:|
| Mandatory ILE61 CD1 mode | 41/45 |
| LEU54 side-chain hydrophobic OR backbone-O HBD modes | 40/45 |
| HIS96 pi modes | 40/45 |
| VAL93 hydrophobic modes | 18/45 |
| LYS94 salt-bridge OR HBA modes | 9/45 |
| All five requirements together | 8/45 |

These are geometric replay counts, not original observation frequencies, an
independent activity panel or full-library retention. The replay is saved at
`data/e054-user-recommendation-rule-audit.json`. Passing poses: 4HBM:0Y7:A:201,
4ERF:0R3:A:201, 4ERE:0R2:A:201, 4QOC:35T:A:201, 4QO4:35S:A:201,
4OCC:2TZ:A:201, 4OBA:2TW:A:501, 4OAS:2SW:A:201.

Inspection found an integration omission: the generic grouped-scoring backend
and edit router supported optional groups, but the recommendation system prompt
still requested only the old fields and permitted an automatically chosen hard
core. The new recommendation is schema-valid under the old policy, but not an
appropriate default for the requested broad-coverage workflow.

The recommendation prompt now explicitly requests optional_groups and explains
their best-member-only scoring and disjoint membership. Automatic recommendations
must leave mandatory_anchors and alternative_groups empty; local validation
enforces this instead of relying only on prompt compliance. The statistical
mandatory_proposal_eligible flag does not authorize the LLM to impose a hard rule.
Explicit human design edits still use the existing support and evidence checks.

Automatic spatial definitions remain empty until a user supplies reviewed source
atom selections and coefficients. Optional feature families are not presented as
complete spatial subpockets. No coordinates, atom selections or nonlinear weights
are invented to fill the missing fields.

Validation: 68 focused tests passed, covering provider workflow, recommendation,
automatic-hard-rule rejection, explicit-human-hard-rule acceptance, optional
families, spatial review requirements and existing consensus behavior.

Workstation continuation: pull the feature branch, restart Chat with the SAME
current workspace and profiles, then submit a NEW `/recommend` in the same
conversation. Existing structure diversity and consensus remain reusable; do not
resume the old recommendation task or create another empty workspace. No source
structures, consensus coordinates or screening coefficients changed in this fix.
