#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:-$(pwd)}
AIDD_PY=${AIDD_PY:-python}
QT9_REPORT=${QT9_REPORT:-$ROOT/data/e020_qt9_validation/validation-report.json}
AGGREGATION_DIR=${E028_AGGREGATION_DIR:-$ROOT/data/e026_8bju_1x8b_aggregation/result}
OUTPUT_DIR=${E028_OUTPUT_DIR:-$ROOT/data/e028_predocking_qc}
RECEPTOR_8BJU=${E028_RECEPTOR_8BJU:-$ROOT/data/e019_query_8bju/8BJU.cif}
RECEPTOR_1X8B=${E028_RECEPTOR_1X8B:-$ROOT/data/e026_query_1x8b/1X8B.cif}
TOP_N=${E028_TOP_N:-100}
CHEMICAL_COMPANION=${CHEMICAL_COMPANION:-}

for path in "$QT9_REPORT" "$AGGREGATION_DIR/manifest.json" \
  "$RECEPTOR_8BJU" "$RECEPTOR_1X8B"; do
  [[ -f "$path" ]] || { echo "Required E028 input not found: $path" >&2; exit 66; }
done

INDEX_CATALOG=$($AIDD_PY -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["result"]["index_catalog"])' \
  "$QT9_REPORT")
ARTIFACT_CATALOG=${ARTIFACT_CATALOG:-$($AIDD_PY -c \
  'import json,sys; print(json.load(open(sys.argv[1]))["artifact_catalog"])' \
  "$INDEX_CATALOG")}
[[ -f "$ARTIFACT_CATALOG" ]] || {
  echo "Artifact catalog not found: $ARTIFACT_CATALOG" >&2; exit 66;
}

CHEMISTRY_ARGS=()
if [[ -n "$CHEMICAL_COMPANION" ]]; then
  [[ -f "$CHEMICAL_COMPANION" ]] || {
    echo "Chemical companion catalog not found: $CHEMICAL_COMPANION" >&2; exit 66;
  }
  CHEMISTRY_ARGS=(--chemical-companion "$CHEMICAL_COMPANION")
fi

$AIDD_PY -m aidd_agent.cli run-predocking-pocket-qc \
  --artifact-catalog "$ARTIFACT_CATALOG" \
  --aggregation-dir "$AGGREGATION_DIR" \
  --output-dir "$OUTPUT_DIR" \
  "${CHEMISTRY_ARGS[@]}" \
  --receptor "8BJU-prepared-v1=$RECEPTOR_8BJU" \
  --receptor "1X8B-prepared-v1=$RECEPTOR_1X8B" \
  --top-n-per-query "$TOP_N" \
  --query-neighborhood 4.0 \
  --close-distance 2.0 \
  --severe-distance 1.5

$AIDD_PY - "$OUTPUT_DIR" "$TOP_N" "$CHEMICAL_COMPANION" <<'PY'
import hashlib
import json
from pathlib import Path
import statistics
import sys

root = Path(sys.argv[1])
top_n = int(sys.argv[2])
chemistry_expected = bool(sys.argv[3])
manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
rows = [json.loads(line) for line in
        (root / "pose-qc.jsonl").read_text(encoding="utf-8").splitlines() if line]
sha256 = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
assert manifest["status"] == "complete"
assert manifest["invariants"]["retrieval_membership_changed"] is False
assert manifest["invariants"]["admission_order_changed"] is False
assert manifest["invariants"]["chemical_topology_available"] is chemistry_expected
assert manifest["invariants"]["point_clouds_are_docking_inputs"] is False
assert len(rows) == manifest["counts"]["selected_tasks"]
assert set(manifest["counts"]["selected_by_query"]) == {
    "8BJU:QT9:A:601", "1X8B:824:A:901"}
assert all(value == top_n for value in manifest["counts"]["selected_by_query"].values())
for output in manifest["outputs"].values():
    assert sha256(Path(output["path"])) == output["sha256"]
assert all(Path(row["point_cloud_pdb"]).is_file() and
           sha256(Path(row["point_cloud_pdb"])) == row["point_cloud_pdb_sha256"]
           for row in rows)
if chemistry_expected:
    assert all(row["vdw_exclusion"] is not None and
               Path(row["chemical_sdf"]).is_file() and
               sha256(Path(row["chemical_sdf"])) == row["chemical_sdf_sha256"]
               for row in rows)
print("E028 accepted coordinate-only execution")
print("selected_by_query=", manifest["counts"]["selected_by_query"])
print("poses_with_close_points=", manifest["counts"]["poses_with_close_points"])
print("poses_with_severe_points=", manifest["counts"]["poses_with_severe_points"])
print("median_centroid_distance_A=", round(statistics.median(
    row["centroid_distance_angstrom"] for row in rows), 3))
print("review=", root / "review.pml")
PY
