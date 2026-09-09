#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:-$(pwd)}
AIDD_PY=${AIDD_PY:-python}
ARTIFACT_CATALOG=${ARTIFACT_CATALOG:-$ROOT/data/e019_artifacts/catalog.json}
CHEMICAL_COMPANION=${CHEMICAL_COMPANION:-/mnt/local/hand/yuzhang/aidd/chemical-companion-v1/catalog.json}
OUTPUT_ROOT=${E031_OUTPUT_DIR:-/mnt/local/hand/yuzhang/aidd/e031-key-interaction-matches-v1}

QT9_QUERY=${E031_QT9_QUERY:-$ROOT/data/e019_query_8bju/gaussian-query-v1.npz}
QT9_RIGID=${E031_QT9_RIGID:-$ROOT/data/e020_qt9_validation/gaussian-staged-v1/refine/merged-scores.npz}
X8B_QUERY=${E031_1X8B_QUERY:-$ROOT/data/e026_query_1x8b/gaussian-query-v1.npz}
X8B_RIGID=${E031_1X8B_RIGID:-$ROOT/data/e026_1x8b_validation/gaussian-staged-v1/refine/merged-scores.npz}

for path in "$ARTIFACT_CATALOG" "$CHEMICAL_COMPANION" \
  "$QT9_QUERY" "$QT9_RIGID" "$X8B_QUERY" "$X8B_RIGID"; do
  [[ -f "$path" ]] || { echo "Required E031 input not found: $path" >&2; exit 66; }
done
mkdir -p "$OUTPUT_ROOT"

run_query() {
  local label=$1 query=$2 rigid=$3 output=$4
  echo "=== E031 $label ==="
  local rigid_hash_before
  rigid_hash_before=$(sha256sum "$rigid" | awk '{print $1}')
  /usr/bin/time -v -o "$OUTPUT_ROOT/$label.time.txt" \
    "$AIDD_PY" -m aidd_agent.cli score-key-interaction-matches \
      --artifact-catalog "$ARTIFACT_CATALOG" \
      --chemical-companion "$CHEMICAL_COMPANION" \
      --query "$query" --rigid-result "$rigid" --output "$output"
  [[ "$(sha256sum "$rigid" | awk '{print $1}')" == "$rigid_hash_before" ]] || {
    echo "Rigid result changed unexpectedly: $rigid" >&2
    exit 70
  }

  "$AIDD_PY" - "$label" "$rigid" "$output" <<'PY'
import json
from pathlib import Path
import sys
import numpy as np

label, rigid_path, match_path = sys.argv[1:]
objective = "atomcentered_anchored_joint"
manifest = json.loads(Path(match_path).with_suffix(".manifest.json").read_text())
with np.load(rigid_path, allow_pickle=False) as rigid, np.load(match_path, allow_pickle=False) as match:
    assert np.array_equal(rigid["global_ids"], match["global_ids"])
    scores = match[f"{objective}__interaction_match_score"]
    baseline = rigid[f"{objective}__objective"]
    order = np.argsort(-scores, kind="stable")
    baseline_order = np.argsort(-baseline, kind="stable")
    baseline_rank = np.empty(len(order), dtype=np.int64)
    baseline_rank[baseline_order] = np.arange(1, len(order) + 1)
    summary = manifest["analysis"][objective]
    print("query:", label)
    print("candidates:", len(scores), "anchors:", manifest["anchors"])
    print("wall_seconds:", manifest["wall_seconds"])
    if label == "8BJU_QT9_A_601":
        reference = 77.99
        overhead_percent = 100.0 * manifest["wall_seconds"] / reference
        print("vs_accepted_77.99s_refine_percent:", overhead_percent)
        print("latency_gate_le_10_percent:", overhead_percent <= 10.0)
    print("interaction_min_median_max:", summary["interaction_min_median_max"])
    print("matched_anchors_median:", summary["matched_anchors_median"])
    print("spearman_vs_gaussian:", summary["spearman_rho_vs_rigid_objective"])
    print("top_k_overlap:", summary["top_k_overlap"])
    print("top_interaction_candidates:")
    for index in order[:10]:
        print({
            "interaction_rank": int(np.where(order == index)[0][0]) + 1,
            "gaussian_rank": int(baseline_rank[index]),
            "global_id": int(match["global_ids"][index]),
            "molecule_id": str(match["molecule_ids"][index]),
            "conformer_id": str(match["conformer_ids"][index]),
            "interaction_match_score": float(scores[index]),
            "anchor_scores": match[f"{objective}__anchor_scores"][index].tolist(),
            "assignments": match[f"{objective}__anchor_assignments"][index].tolist(),
        })
PY
}

run_query "8BJU_QT9_A_601" "$QT9_QUERY" "$QT9_RIGID" \
  "$OUTPUT_ROOT/8BJU_QT9_A_601.interaction-matches.npz"
run_query "1X8B_824_A_901" "$X8B_QUERY" "$X8B_RIGID" \
  "$OUTPUT_ROOT/1X8B_824_A_901.interaction-matches.npz"

echo "E031 outputs: $OUTPUT_ROOT"
