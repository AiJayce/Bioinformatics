#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

mkdir -p data/raw data/models data/reference data/exports

need() {
  local name="$1"
  local value="${!name:-}"
  if [[ -z "$value" ]]; then
    echo "[ERROR] Missing required URL env var: $name"
    exit 1
  fi
}

download_file() {
  local url="$1"
  local out="$2"
  mkdir -p "$(dirname "$out")"
  echo "[INFO] Download $url -> $out"
  curl -L --fail --retry 3 "$url" -o "$out"
}

extract_archive() {
  local archive="$1"
  local target_dir="$2"
  mkdir -p "$target_dir"
  case "$archive" in
    *.tar.gz|*.tgz) tar -xzf "$archive" -C "$target_dir" ;;
    *.zip) unzip -o "$archive" -d "$target_dir" ;;
    *)
      echo "[ERROR] Unsupported archive format: $archive"
      exit 1
      ;;
  esac
}

# Required 3 inputs (user-requested core set)
need RAW_MERSCOPE_URL
need RAW_SCRNA_URL
need VAE_MODEL_ARCHIVE_URL

download_file "$RAW_MERSCOPE_URL" data/raw/merscope_raw.h5ad
download_file "$RAW_SCRNA_URL" data/raw/scrna_reference_raw.h5ad

VAE_ARCHIVE="data/models/vae_model_bundle.tmp"
download_file "$VAE_MODEL_ARCHIVE_URL" "$VAE_ARCHIVE"
extract_archive "$VAE_ARCHIVE" data/models
rm -f "$VAE_ARCHIVE"

# Optional but recommended pipeline dependencies
if [[ -n "${PLAQUES_URL:-}" ]]; then
  download_file "$PLAQUES_URL" data/raw/Plaques_Object.h5ad
fi

if [[ -n "${ATLAS_ARCHIVE_URL:-}" ]]; then
  ATLAS_ARCHIVE="data/reference/atlas_bundle.tmp"
  download_file "$ATLAS_ARCHIVE_URL" "$ATLAS_ARCHIVE"
  extract_archive "$ATLAS_ARCHIVE" data/reference
  rm -f "$ATLAS_ARCHIVE"
fi

if [[ -n "${MERSCOPE_REGIONS_ARCHIVE_URL:-}" ]]; then
  REGION_ARCHIVE="data/raw/merscope_regions_bundle.tmp"
  download_file "$MERSCOPE_REGIONS_ARCHIVE_URL" "$REGION_ARCHIVE"
  extract_archive "$REGION_ARCHIVE" data/raw
  rm -f "$REGION_ARCHIVE"
fi

if [[ -n "${EXPORTS_ARCHIVE_URL:-}" ]]; then
  EXPORTS_ARCHIVE="data/exports/exports_bundle.tmp"
  download_file "$EXPORTS_ARCHIVE_URL" "$EXPORTS_ARCHIVE"
  extract_archive "$EXPORTS_ARCHIVE" data/exports
  rm -f "$EXPORTS_ARCHIVE"
fi

if ! find data/models -type f -name 'vae_model.pt' | grep -q .; then
  echo "[ERROR] VAE model import failed: no vae_model.pt found under data/models"
  exit 1
fi

echo "[INFO] Data fetch complete."
