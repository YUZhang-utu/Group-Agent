# Contact-first search from Chat

The new `budget` workflow runs molecule-level ANN retrieval, expands all stored
conformers for the retrieved molecules, computes physical pocket exclusion and
same-pose contact/Gaussian scores, and exports the first molecule page.

Within each template, both pose selection and molecule ranking use contact score
first and Gaussian score only for ties. The old 0.7/0.3 composite is diagnostic.
Template ranks are fused by RRF. Any-positive-contact molecules precede the
zero-contact reserve; within each tier, per-template quota membership precedes
the reserve and RRF determines order. The complete scored ranking is saved.

Defaults: 1,000,000 retrieved unique molecules, 100,000 exported unique molecules,
100,000 per-template quota, RRF k=60. These are computational budgets, not
validated activity cutoffs. If fewer molecules have a physical surviving pose,
export reports the shortfall rather than duplicating molecules. Retrieval uses
the installed USRCAT FAISS index, not a newly validated pharmacophore index.
The current 8,318,351-molecule library is not a 100-million-molecule library;
1% would not even supply the requested first 100,000 molecules.

## Update the workstation

After current jobs finish, update the existing checkout and restart the Chat
server using the SAME storage root, runtime and provider configuration as before:

```bash
conda activate aidd-workstation
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin feature/structure-guided-chat
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export AIDD_ASSIGNMENT_BACKEND=numba
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -u -m aidd_agent.chat_web \
  --storage-root /mnt/local/hand/yuzhang/aidd/e054-chat-20260923-104806 \
  --runtime /mnt/local/hand/yuzhang/aidd/config/prompt-runtime-apptainer.local.json \
  --allow-compute --port 8765
```

Preserve an existing `--llm-config-dir` argument or `AIDD_LLM_CONFIG_DIR` if used.
The trusted runtime needs `search.batch`, `search.workers`, and
`search.refine_chunk`. Existing runtime search fields remain supported.
Use the existing MDM2 conversation containing a completed consensus recommendation
or design. Natural-language request:

> Run budget search for this MDM2 project: retrieve 1000000 unique molecules,
> rank contacts first with shape only for ties, fuse template ranks, and export
> the first 100000 molecules with names and original and posed MOL2 files.

Deterministic equivalent: `/budget`, optionally followed by the source Chat task
ID. `/status` reads progress and `/results` returns the completed page paths.
After completion, request the next 100000 molecules or use `/budget_page`.
Default pagination starts immediately after the last completed delivery. Explicit
source task IDs scope pagination to that delivery; retain the latest source to
avoid returning to an older page. Resume an interrupted task with `/resume`.

Each delivery contains `molecules.csv`, `molecules.jsonl`, `original/`, `posed/`
and an integrity report. Original MOL2 names, hydrogens and bonds are preserved;
posed files change coordinates only. These are not protonation-prepared docking
inputs. Export stops for source identity/hash mismatches and reports chemistry
review flags. Later pages reuse the ranking without repeating pose scoring.

Existing water/metal-mediated anchor roles are removed in a recorded protein-only
copy. Original recommendations/designs remain unchanged. Source coordinates,
target identity and templates come from the owned project. Reviewed occupancy
definitions are optional and must match target and coordinate frame. The MDM2
fixture is never reused for another target; absent definitions produce unknown
occupancy. New targets still need their own completed consensus preparation.

Local regression does not establish live provider routing, workstation throughput,
MDM2 ANN recall or activity enrichment. Start the workstation job to obtain those
execution results; completion today cannot be inferred from the old small pilot.
Old E059 runs cannot resume under changed ranking/code fingerprints; start a new
budget task. Do not change code or runtime during that sealed task.
