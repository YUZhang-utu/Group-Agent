# E144: Explicit PDB inspection

## Protocol

Hypothesis: an allowlisted, target-independent PDB inspection action can report
deposited evidence without requiring prior isoform or receptor selection.

Implement `pdb_inspect` with a single PDB ID. Preserve the downloaded mmCIF and
its digest, save a JSON report in the task directory, and expose it through the
existing project task/artifact tools. Do not run preparation or inference.

Validate with synthetic mmCIF cases covering two asymmetric-unit protein chains
versus separate monomer assembly annotations, sequence/coordinate differences,
missing metadata, author/label identifiers, ions, modified polymer residues,
and database mapping/differences. Verify plan validation and executor integration.
Use a public 4TQ9 download for a separate live acquisition/parser check if network
access permits. This does not establish live-provider routing acceptance.

## Acceptance

Missing metadata remains explicitly unreported. Database-associated accession
is not an independently verified isoform. Assembly annotations are not inferred
from asymmetric-unit chain counts. No ligand is automatically selected for search.
Existing task reports and project ownership remain authoritative persistence.

## Results

Confirmatory local validation: 51 tests passed across inspection, general prompt,
project guidance and domain-agent suites. A second 60-test group had 59 passes
and one obsolete session-isolation assertion. Updated that assertion to verify
same-project sharing and different-project isolation; all three tests in that
module passed on rerun. Total unique relevant cases passing: 111.

English-content and whitespace checks passed. Live 4TQ9 acquisition was attempted
but not completed: initial test-directory access and then HTTPS socket access
were denied by the environment; escalated attempts exited without diagnostics.
No live acquisition/parser or cloud-provider routing acceptance is claimed.

The general clarification-to-revised-plan lifecycle remains a separate open issue.
