# Data Layout (Repository-local)

This repository is configured to read all inputs from `data/`.

Required core files:
- `data/raw/merscope_raw.h5ad`
- `data/raw/scrna_reference_raw.h5ad`
- `data/models/**/vae_model.pt` (plus associated VAE outputs such as latent/ridge files)

Recommended additional files for full figure parity:
- `data/raw/Plaques_Object.h5ad`
- `data/raw/merscope_regions/.../region_*.h5ad`
- `data/reference/DAPI_downsampling/...`
- `data/exports/DESeq2/...`
- `data/exports/Nichnet/...`
- `data/exports/fGSEA/...`

Use `./fetch_required_data.sh` to download these files from your own URLs/releases.
