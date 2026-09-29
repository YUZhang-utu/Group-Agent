# Boundary candidates and continuous property/steric review

Use broad continuous side-chain properties, not exact sequence equality. Steric
size must remain a separate signal: equal hydrophobic annotation does not imply
equal local shape, reach or proximal branching. This implementation does not
yet provide a calibrated final property-based partition or production dispatcher.

## Run the whole-boundary inspection now

The completed backbone build is sufficient. No new MOL2 extraction is required.
Update feature/structure-guided-chat, then run from the repository root:

```bash
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
BUILD=/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/backbone
python -m aidd_agent.boundary_block_candidates \
  --build "$BUILD" \
  --diagnostic "${BUILD}-fragmentation-review" \
  --output "${BUILD}-boundary-candidates-v1"
```

This inspects all 1,384,339 boundary descriptors by indexed lookups. It writes
one row per conformer to candidates.sqlite, including original class, candidate
definite classes, allowed directed rotations and reasons for unresolved cases.
It can take appreciable I/O time; progress is printed every 10,000 records.
It does not rewrite existing membership. Use a fresh output if interrupted.

report.json contains the actual 5-degree angle histogram and candidate coverage
at deviations 35/45/60/75 degrees from cis/trans centres. These are sensitivity
settings, not adopted thresholds. The reported regular coverage is an upper
bound before property and geometry checks, not achieved final coverage.
No cutoff is selected just to exceed 95%.

The old descriptor does not expose chirality directly. The tool verifies finite
R/S/achiral sequences against the original class hash, with a bounded trial
budget. Unresolved or unsupported signatures remain special. Known cis/trans
positions cannot change; only directed cyclic rotations are allowed. The
90-degree midpoint is not assigned. Candidate classes must already exist as
definite classes; no union of definite classes through boundary records occurs.

## Side-chain properties and steric review

macrocycle_property_profiles provides continuous property and steric diagnostics:

- Existing chemistry features: side-chain heavy-atom size proxy, formal charge,
  cyclic N/N-methyl features and donor/acceptor/aromatic/hydrophobic counts.
- Existing typed features: local feature-coordinate first/second moments.
- New sample-level steric features: actual all-side-heavy-atom VDW-expanded
  local-axis extents, reach, covariance eigenvalues and proximal branch degree.

These are coordinate-dependent steric proxies, not exact union volumes, clash
energies or complete torsional barriers. Properties and sterics are compared
separately; no arbitrary combined weight or hard exact-property class is imposed.
The legacy typed and chemistry builds remain unchanged and reusable.

For an already exported review MOL2, compute at most 200 profiles:

```bash
python -m aidd_agent.macrocycle_property_profiles \
  --mol2 REVIEW_SAMPLE.mol2 --output REVIEW_PROPERTIES.jsonl --limit 200
```

This command is for selected review molecules, not the entire raw MOL2 library.
The output must not already exist. A failed strict molecule parse stops the
review; it is not silently accepted. Pairwise property_difference compares
corresponding residues under a verified rotation; candidate identity/chirality
admission must precede this comparison.

Report sibling receipt existence automatically, or check on the workstation:

```bash
for v in chemistry typed; do
  test -f "${BUILD%/backbone}/$v/report.json" && echo "$v receipt exists"
  test -f "${BUILD%/backbone}/$v/descriptors.sqlite" && echo "$v descriptors exist"
done
```

File existence is not completion or provenance validation. Inspect status and
source identity before reuse. Do not rerun all variants merely to answer this.

## Wrapper correction

The previous run_e080_macrocycle_blocks.sh used set -e, while successful
structural validation with review_required returns 2. This could stop the
variant loop after backbone. A new explicit --allow-review-for-offline option
permits continuation only for complete, error-free structural validation;
release_gate remains review_required. Default validator exit behaviour remains
unchanged. Actual structural failures still stop the wrapper.

## Local evidence and outstanding work

Small existing backbone fixture: 24 conformers, 3 boundary records; all 9 group
chirality signatures recovered. Two boundary records have angular candidates
at 35 degrees, all three at 45. These are local fixture results only.
Unit tests cover angle/state consistency, chirality and cis/trans preservation,
midpoint rejection, cyclic rotation, rigid-transform invariance, side-chain
coordinate response and offline validator exit policy.

Workstation angle counts and chemistry/typed availability remain pending.
Final property-compatible assignment, special-set search timing and held-out
recall still require these results and evaluation. Do not interpret this
candidate tool or the earlier queue plan as a finished scientific partition.

Reference: RDKit documentation describes property and 3D shape descriptors:
https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-3d-descriptors
The local steric profile here is explicitly defined above, not an RDKit union-volume call.
