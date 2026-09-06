#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 7 || $# -gt 8 ]]; then
  echo "Usage: bash scripts/validate_qt9_pharmacophore.sh ARTIFACT_CATALOG FAISS_INDEX_DIR QUERY_MANIFEST MMCIF CCD PHARMACOPHORE_INDEX_DIR OUTPUT_DIR [CCD_ID]" >&2
  exit 64
fi

ARTIFACT_CATALOG=$1
FAISS_INDEX_DIR=$2
QUERY_MANIFEST=$3
MMCIF=$4
CCD=$5
PHARMACOPHORE_INDEX_DIR=$6
OUTPUT_DIR=$7
CCD_ID=${8:-QT9}

for required_file in "$ARTIFACT_CATALOG" "$QUERY_MANIFEST" "$MMCIF" "$CCD"; do
  if [[ ! -f "$required_file" ]]; then
    echo "Required file not found: $required_file" >&2
    exit 66
  fi
done
for required_file in "$FAISS_INDEX_DIR/manifest.json" "$FAISS_INDEX_DIR/transform.npz"; do
  if [[ ! -f "$required_file" ]]; then
    echo "Required FAISS file not found: $required_file" >&2
    exit 66
  fi
done

mkdir -p "$PHARMACOPHORE_INDEX_DIR" "$OUTPUT_DIR"

python -m aidd_agent.cli validate-pharmacophore-retrieval \
  --artifact-catalog "$ARTIFACT_CATALOG" \
  --pharmacophore-index-dir "$PHARMACOPHORE_INDEX_DIR" \
  --query-manifest "$QUERY_MANIFEST" \
  --output-dir "$OUTPUT_DIR" \
  --faiss-index-dir "$FAISS_INDEX_DIR" \
  --mmcif "$MMCIF" \
  --ccd "$CCD" \
  --ccd-id "$CCD_ID" \
  --faiss-k 100000 \
  --nprobe 256 \
  --bin-width 0.5 \
  --max-distance 20.0 \
  2>&1 | tee "$OUTPUT_DIR/validation.log"

python - "$OUTPUT_DIR/validation-report.json" <<'PY'
import json
import sys
from pathlib import Path

path = Path(sys.argv[1])
report = json.loads(path.read_text(encoding="utf-8"))
print("\n=== QT9 local pharmacophore validation summary ===")
print("accepted:", report["accepted"])
print("library conformers:", report["index"]["conformers"])
print("indexed shards:", report["index"]["shards"])
print("supported anchors:", report["query"]["supported_anchors"])
print("query pairs:", report["query"]["pairs"])
print("FAISS returned:", report["faiss"]["returned"])
print("FAISS query seconds:", report["faiss"]["query_seconds"])
print("tier counts:", report["validation"]["counts"])
print("timings:", report["timings"])
failed = [name for name, passed in report["validation"]["checks"].items() if not passed]
print("failed checks:", failed or "none")
if not report["accepted"]:
    raise SystemExit(2)
PY
