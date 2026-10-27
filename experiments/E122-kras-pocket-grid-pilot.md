# E122: rapid protein-grid pilot within retained KRAS pools

Protocol written before execution. Use existing sealed 7RPZ and 12FF receptor
coordinates, including the same-chain GDP/Mg already in the geometry package.
Do not replace G12V with 9ZPA. Keep original typed same-pose contact rules.

Sample 10,000 unique mapped pool IDs per target uniformly without replacement
using a fixed RNG seed. Sample original pool membership, not top scoring IDs.
Read only compact molecule/conformer identity arrays once to map sampled IDs;
do not traverse full-library coordinates. Evaluate up to three existing conformers
per sampled ID, taking the first by global ID, and the original bounded seed set.
Share reads across targets. Missing mappings are unresolved and reported.

Build a 0.75 A nearest-receptor-atom distance grid in the reference-located box
with 8 A padding. The ligand locates the site only: pocket burial is derived
from protein occupancy along six Cartesian rays (2,4,6,8 A), not ligand overlap.
At least three obstructed rays defines the exploratory buried-region proxy.
This is not a complete cavity-volume calculation or an atom-type force field.
Allow solvent-exposed atoms; no requirement that every atom be buried.

For direct-contact poses, interpolate grid distances and burial flags. Rank by
burial fraction minus mean squared overlap below 2.2 A minus twice the fraction
below 1.2 A. Recheck the best eight poses per conformer using exact receptor
distances and exact six-ray tests. Report three predeclared exploratory profiles
on exact metrics: (max severe fraction, max mean squared overlap, min burial)
lenient=(0.10,0.60,0.15), balanced=(0.05,0.40,0.25), strict=(0,0.25,0.35).
Profiles permit some geometric conflict and are prioritization, not valid docking
poses or proof of impossibility. Keep misses archived; do not delete source IDs.

Controls: native direct contacts, native exact metrics for all profiles, and a
ligand translated into a receptor atom must fail the severe-distance check.
Report rather than hide any profile failing native recovery. Main sample proceeds
only if both native controls retain the lenient profile and original contacts.
Compare interpolated distances with exact distances and expose discrepancy.

Outputs: sampled IDs with conformers, per-molecule best pose and profile flags,
target retention estimates with binomial Wilson intervals, compute time separate
from setup/index time. No extrapolated hit or enrichment claim. No automatic
million-ID truncation or production export in this pilot. Use empirical retention
and throughput to choose the next run. Controls do not establish candidate recall.
