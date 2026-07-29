#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

usage() {
  cat <<'EOF'
Usage:
  ./run_e2e.sh [cores]

Optional environment variables:
  SNAKEMAKE_BIN   Path to snakemake binary (default: snakemake)
  CONFIG_FILE     Path to config file (default: config/config.run.yaml)
  CORES           Number of cores (default: first arg or 4)
  DRY_RUN         If set to 1, run with -n only
  FETCH_DATA      If set to 1, run ./fetch_required_data.sh before Snakemake
  SCVI_PYTHON     Override runtime.scvi_python in resolved config
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

SNAKEMAKE_BIN="${SNAKEMAKE_BIN:-snakemake}"
CONFIG_FILE="${CONFIG_FILE:-config/config.run.yaml}"
CORES="${CORES:-${1:-4}}"

if [[ ! -f "$CONFIG_FILE" ]]; then
  cp config/config.portable.yaml "$CONFIG_FILE"
  echo "[INFO] Created $CONFIG_FILE from config/config.portable.yaml"
fi

RESOLVED_CONFIG="config/config.resolved.yaml"
cp "$CONFIG_FILE" "$RESOLVED_CONFIG"

# Backward-compatible substitution when ALZ_DATA_ROOT is provided.
if [[ -n "${ALZ_DATA_ROOT:-}" ]]; then
  sed -i "s#\${ALZ_DATA_ROOT}#${ALZ_DATA_ROOT}#g" "$RESOLVED_CONFIG"
fi

if [[ -n "${SCVI_PYTHON:-}" ]]; then
  if grep -q '^runtime:' "$RESOLVED_CONFIG"; then
    sed -i "s#^  scvi_python:.*#  scvi_python: ${SCVI_PYTHON}#" "$RESOLVED_CONFIG"
  else
    {
      echo ""
      echo "runtime:"
      echo "  scvi_python: ${SCVI_PYTHON}"
    } >> "$RESOLVED_CONFIG"
  fi
fi

SCVI_PY="$(awk '/^runtime:/{flag=1;next}/^[^ ]/{flag=0}flag && $1=="scvi_python:"{print $2}' "$RESOLVED_CONFIG" | tail -n 1)"
if [[ -z "$SCVI_PY" ]]; then
  SCVI_PY="python"
fi

"$SCVI_PY" - <<'PY'
import importlib.util, sys
ok = importlib.util.find_spec('scvi') is not None
if not ok:
    raise SystemExit("[ERROR] scvi-tools is not installed in runtime.scvi_python")
print("[INFO] scvi-tools check: OK")
PY

if [[ "${FETCH_DATA:-0}" == "1" ]]; then
  if [[ ! -x "./fetch_required_data.sh" ]]; then
    echo "[ERROR] fetch_required_data.sh not found or not executable"
    exit 1
  fi
  ./fetch_required_data.sh
fi

if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "[INFO] Dry-run mode"
  exec "$SNAKEMAKE_BIN" -s Snakefile --configfile "$RESOLVED_CONFIG" -n -p
fi

exec "$SNAKEMAKE_BIN" -s Snakefile --configfile "$RESOLVED_CONFIG" --cores "$CORES" -p
