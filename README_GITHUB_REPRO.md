# Git Clone Reproducibility Guide

This pipeline is configured to run from repository-local `data/` paths only.
No server-specific absolute paths are required.

## One-command end-to-end run

After clone, fetch data to `data/` and run:

```bash
FETCH_DATA=1 ./run_e2e.sh 4
```

Dry-run check:

```bash
DRY_RUN=1 ./run_e2e.sh
```

## 1) Clone

```bash
git clone <YOUR_REPO_URL>
cd Snakemake
```

## 2) Prepare environment

Use the same Snakemake Python environment you used before, or create a new conda env with required packages.

```bash
conda create -n ccf_snakemake python=3.11 -y
conda activate ccf_snakemake
pip install snakemake scanpy anndata pandas numpy matplotlib seaborn pynrrd scipy matplotlib-venn
```

R dependency for fgsea step:

```bash
# inside R
# install.packages(c("BiocManager", "ggplot2", "stringr"))
# BiocManager::install(c("fgsea", "msigdbr"))
```

## 3) Put data/models under one root directory

Repository-local data path is used:

```bash
data/
```

To auto-download core assets, define your own URLs (for example GitHub Release assets, object storage, or institutional file server):

```bash
export RAW_MERSCOPE_URL="https://.../merscope_raw.h5ad"
export RAW_SCRNA_URL="https://.../scrna_reference_raw.h5ad"
export VAE_MODEL_ARCHIVE_URL="https://.../vae_model_bundle.tar.gz"

# optional for full notebook parity
export PLAQUES_URL="https://.../Plaques_Object.h5ad"
export ATLAS_ARCHIVE_URL="https://.../DAPI_downsampling.tar.gz"
export MERSCOPE_REGIONS_ARCHIVE_URL="https://.../merscope_regions.tar.gz"
export EXPORTS_ARCHIVE_URL="https://.../exports_bundle.tar.gz"
```

Then run:

```bash
chmod +x fetch_required_data.sh run_e2e.sh
FETCH_DATA=1 ./run_e2e.sh 4
```

## 4) Config

Default config already points to repository-local paths.

- Main config: `config/config.yaml`
- Runner template copy: `config/config.portable.yaml`

## 5) Run

Dry-run:

```bash
snakemake -s Snakefile --configfile config/config.run.yaml -n -p
```

Actual run:

```bash
snakemake -s Snakefile --configfile config/config.run.yaml --cores 4 -p
```

Equivalent one-command wrapper:

```bash
./run_e2e.sh 4
```

## scVI and VAE reproducibility scope

- This workflow runs **raw query -> scVI/SCANVI annotation** in pipeline rule `scvi_annotate_raw`.
- Provide both:
  - `source.raw_merscope_h5ad` (query raw data)
  - `source.scvi_reference_h5ad` (labeled reference data with `scvi.label_key`)
- VAE is **not retrained**. The workflow imports pretrained artifacts from `VAE_model/` (including `vae_model.pt`, latent/ridge outputs) for downstream plots.

If `scvi-tools` is not in your active Python, set `runtime.scvi_python` in the config to the exact interpreter path.
