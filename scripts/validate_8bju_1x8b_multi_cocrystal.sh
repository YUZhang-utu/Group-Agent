#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:-$(pwd)}
AIDD_PY=${AIDD_PY:-python}
QT9_VALIDATION=${QT9_VALIDATION:-$ROOT/data/e020_qt9_validation}
QT9_RUN=${QT9_RUN:-$QT9_VALIDATION/gaussian-staged-v1}
QUERY_DIR=${QUERY_DIR:-$ROOT/data/e026_query_1x8b}
VALIDATION_DIR=${VALIDATION_DIR:-$ROOT/data/e026_1x8b_validation}
AGGREGATION_DIR=${AGGREGATION_DIR:-$ROOT/data/e026_8bju_1x8b_aggregation}

QT9_REPORT=$QT9_VALIDATION/validation-report.json
QT9_RESULT=$QT9_RUN/refine/merged-scores.npz
QT9_RESULT_MANIFEST=$QT9_RUN/refine/merged-scores.manifest.json
for path in "$QT9_REPORT" "$QT9_RESULT" "$QT9_RESULT_MANIFEST"; do
  [[ -f "$path" ]] || { echo "Required QT9 file not found: $path" >&2; exit 66; }
done

INDEX_CATALOG=${INDEX_CATALOG:-$($AIDD_PY -c 'import json,sys; print(json.load(open(sys.argv[1]))["result"]["index_catalog"])' "$QT9_REPORT")}
ARTIFACT_CATALOG=${ARTIFACT_CATALOG:-$($AIDD_PY -c 'import json,sys; print(json.load(open(sys.argv[1]))["artifact_catalog"])' "$INDEX_CATALOG")}
FAISS_INDEX_FILE=${FAISS_INDEX_FILE:-$($AIDD_PY -c 'import json,sys; print(json.load(open(sys.argv[1]))["faiss"]["index"])' "$QT9_REPORT")}
QT9_QUERY_ID=$($AIDD_PY -c 'import json,sys; print(json.load(open(sys.argv[1]))["query"]["query_id"])' "$QT9_REPORT")
FAISS_INDEX_DIR=$(dirname "$FAISS_INDEX_FILE")
PHARMA_INDEX_DIR=$(dirname "$INDEX_CATALOG")
WORKERS=${WORKERS:-16}

mkdir -p "$QUERY_DIR" "$VALIDATION_DIR" "$AGGREGATION_DIR"
MMCIF=$QUERY_DIR/1X8B.cif
CCD=$QUERY_DIR/824.cif
QUERY_MANIFEST=$QUERY_DIR/query_manifest.json
GAUSSIAN_QUERY=$QUERY_DIR/gaussian-query-v1.npz
GAUSSIAN_RUN=$VALIDATION_DIR/gaussian-staged-v1
MULTI_PLAN=$AGGREGATION_DIR/multi-cocrystal-plan.json
MULTI_OUTPUT=$AGGREGATION_DIR/result

if [[ ! -f "$MMCIF" ]]; then
  curl --fail --location --retry 3 --output "$MMCIF.partial" \
    https://files.rcsb.org/download/1X8B.cif
  mv "$MMCIF.partial" "$MMCIF"
fi
if [[ ! -f "$CCD" ]]; then
  curl --fail --location --retry 3 --output "$CCD.partial" \
    https://files.rcsb.org/ligands/download/824.cif
  mv "$CCD.partial" "$CCD"
fi

$AIDD_PY - "$MMCIF" <<'PY'
from pathlib import Path
import sys
from aidd_agent.chemistry_prep import enumerate_ligand_instances
rows = enumerate_ligand_instances(Path(sys.argv[1]), ["824"])
assert len(rows) == 1, rows
row = rows[0]
assert (row["chain_id"], row["residue_number"]) == ("A", "901"), row
print("verified ligand: 1X8B:824:A:901; atoms=", len(row["atoms"]))
PY

if [[ ! -f "$QUERY_MANIFEST" ]]; then
  $AIDD_PY -m aidd_agent.cli extract-query-anchors \
    --mmcif "$MMCIF" --ccd "$CCD" --ccd-id 824 \
    --query-id 1X8B:824:A:901 --output "$QUERY_MANIFEST"
else
  echo "Reusing query manifest: $QUERY_MANIFEST"
fi

RETRIEVAL_READY=false
if [[ -f "$VALIDATION_DIR/validation-report.json" \
      && -f "$VALIDATION_DIR/pharmacophore-hits.npz" \
      && -f "$VALIDATION_DIR/pharmacophore-hits.manifest.json" ]]; then
  if $AIDD_PY -c 'import json,sys; raise SystemExit(0 if json.load(open(sys.argv[1])).get("accepted") is True else 1)' \
      "$VALIDATION_DIR/validation-report.json"; then
    RETRIEVAL_READY=true
  fi
fi
if [[ "$RETRIEVAL_READY" == true ]]; then
  echo "Reusing accepted retrieval: $VALIDATION_DIR/pharmacophore-hits.npz"
else
  $AIDD_PY -m aidd_agent.cli validate-pharmacophore-retrieval \
    --artifact-catalog "$ARTIFACT_CATALOG" \
    --pharmacophore-index-dir "$PHARMA_INDEX_DIR" \
    --query-manifest "$QUERY_MANIFEST" \
    --output-dir "$VALIDATION_DIR" \
    --faiss-index-dir "$FAISS_INDEX_DIR" \
    --mmcif "$MMCIF" --ccd "$CCD" --ccd-id 824 \
    --faiss-k 100000 --nprobe 256 --bin-width 0.5 --max-distance 20.0
fi

if [[ ! -f "$GAUSSIAN_QUERY" ]]; then
  $AIDD_PY -m aidd_agent.cli prepare-gaussian-query \
    --mmcif "$MMCIF" --ccd "$CCD" --query-manifest "$QUERY_MANIFEST" \
    --output "$GAUSSIAN_QUERY"
else
  echo "Reusing Gaussian query: $GAUSSIAN_QUERY"
fi

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
$AIDD_PY -m aidd_agent.cli run-scaled-gaussian-reranking \
  --artifact-catalog "$ARTIFACT_CATALOG" \
  --query "$GAUSSIAN_QUERY" \
  --candidate-schedule "$VALIDATION_DIR/pharmacophore-hits.npz" \
  --output-dir "$GAUSSIAN_RUN" --stage all --workers "$WORKERS" \
  --coarse-chunk-size 2000 --refine-chunk-size 250 \
  --top-n-per-objective 5000 --max-pair-seeds 512

$AIDD_PY - "$MULTI_PLAN" "$QT9_RESULT" "$QT9_RESULT_MANIFEST" \
  "$GAUSSIAN_RUN/refine/merged-scores.npz" \
  "$GAUSSIAN_RUN/refine/merged-scores.manifest.json" "$QT9_QUERY_ID" <<'PY'
import json
from pathlib import Path
import sys
plan = {
    "format": "aidd-multi-cocrystal-query-plan", "version": 1,
    "library_id": "LIB-AFA68EE6888C", "queries": [
        {"query_id": sys.argv[6], "receptor_id": "8BJU-prepared-v1",
         "site_id": "WEE1-ATP-site", "result": sys.argv[2],
         "result_manifest": sys.argv[3],
         "metadata": {"pdb_id": "8BJU", "ccd_id": "QT9"}},
        {"query_id": "1X8B:824:A:901", "receptor_id": "1X8B-prepared-v1",
         "site_id": "WEE1-ATP-site", "result": sys.argv[4],
         "result_manifest": sys.argv[5],
         "metadata": {"pdb_id": "1X8B", "ccd_id": "824"}},
    ]}
path = Path(sys.argv[1]); path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(json.dumps(plan, indent=2), encoding="utf-8")
print("wrote", path)
PY

$AIDD_PY -m aidd_agent.cli aggregate-multi-cocrystal-search \
  --plan "$MULTI_PLAN" --output-dir "$MULTI_OUTPUT" \
  --primary-objective atomcentered_anchored_joint \
  --top-conformers-per-query 3 --per-query-quota 1000 \
  --consensus-quota 5000 --global-limit-per-site 20000 --rrf-k 60

$AIDD_PY - "$VALIDATION_DIR/validation-report.json" "$GAUSSIAN_RUN/run-manifest.json" \
  "$MULTI_OUTPUT" "$MMCIF" "$CCD" "$QUERY_MANIFEST" "$MULTI_PLAN" \
  "$QT9_QUERY_ID" "$AGGREGATION_DIR/e026-validation-summary.json" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

retrieval = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
gaussian = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
root = Path(sys.argv[3])
manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
load = lambda name: [json.loads(line) for line in
                     (root / name).read_text(encoding="utf-8").splitlines() if line]
evidence = load("query-molecule-evidence.jsonl")
summaries = load("molecule-summary.jsonl")
admissions = load("docking-admission.jsonl")
tasks = load("docking-tasks.jsonl")

assert retrieval["accepted"] is True
assert retrieval["validation"]["checks"]["external_l1_preserved"] is True
assert gaussian["status"] == "complete"
assert gaussian["stages"]["coarse"]["candidates"] == retrieval["validation"]["counts"]["union"]
assert gaussian["stages"]["refine"]["candidates"] > 0
assert manifest["status"] == "complete"
assert manifest["counts"]["queries"] == 2
assert manifest["counts"]["sites"] == 1
inv = manifest["invariants"]
assert inv["cross_query_policy"] == "site-scoped union"
assert inv["raw_scores_averaged"] is False
assert inv["query_protected_quotas"] is True
assert inv["unique_admission_keys"] is True
assert inv["all_tasks_reference_admission"] is True

evidence_union = {(row["site_id"], row["molecule_id"]) for row in evidence}
summary_keys = {(row["site_id"], row["molecule_id"]) for row in summaries}
assert evidence_union == summary_keys
support = {(row["site_id"], row["molecule_id"]): set() for row in evidence}
for row in evidence:
    support[(row["site_id"], row["molecule_id"])].add(row["query_id"])
task_support = {}
for row in tasks:
    task_support.setdefault((row["site_id"], row["molecule_id"]), set()).add(row["query_id"])
for row in admissions:
    key = (row["site_id"], row["molecule_id"])
    assert task_support[key] == support[key]

shared = sum(len(value) == 2 for value in support.values())
q1_only = sum(value == {sys.argv[8]} for value in support.values())
q2_only = sum(value == {"1X8B:824:A:901"} for value in support.values())
summary = {
    "format": "aidd-e026-real-dual-cocrystal-validation", "version": 1,
    "accepted": True, "retrieval_1x8b_accepted": True,
    "queries": 2, "sites": 1, "union_molecules": len(support),
    "shared_molecules": shared, "8bju_only_molecules": q1_only,
    "1x8b_only_molecules": q2_only, "admitted_molecules": len(admissions),
    "docking_tasks": len(tasks), "raw_scores_averaged": False,
    "1x8b_supported_anchors": retrieval["query"]["supported_anchors"],
    "1x8b_anchor_pairs": retrieval["query"]["pairs"],
    "1x8b_retrieval_candidates": retrieval["validation"]["counts"]["union"],
    "1x8b_refined_conformers": gaussian["stages"]["refine"]["candidates"],
    "input_sha256": {
        Path(path).name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
        for path in sys.argv[4:8]
    },
}
Path(sys.argv[9]).write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
PY

echo "E026 complete: $AGGREGATION_DIR/e026-validation-summary.json"
