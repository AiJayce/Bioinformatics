import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import yaml


def ensure_parent(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)


def save_placeholder(path, title):
    ensure_parent(path)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, title, ha='center', va='center')
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


src_h5ad = snakemake.input.source
out_h5ad = snakemake.output.h5ad
qc_png = snakemake.output.qc
umap_png = snakemake.output.umap
leiden_png = snakemake.output.leiden

cfg = getattr(snakemake, 'config', None)
if not isinstance(cfg, dict) or len(cfg) == 0:
    with open(snakemake.input.config, 'r', encoding='utf-8') as f:
        cfg = yaml.safe_load(f)

pp_cfg = cfg.get('scrna_preprocess', {})

min_genes_per_cell = pp_cfg.get('min_genes_per_cell', 200)
max_genes_per_cell = pp_cfg.get('max_genes_per_cell', None)
max_pct_counts_mt = pp_cfg.get('max_pct_counts_mt', 20.0)
min_counts_per_cell = pp_cfg.get('min_counts_per_cell', None)
min_cells_per_gene = pp_cfg.get('min_cells_per_gene', 3)
target_sum = pp_cfg.get('target_sum', 10000)
n_top_genes = pp_cfg.get('n_top_genes', 2000)
hvg_flavor = pp_cfg.get('hvg_flavor', 'seurat')
n_pcs = pp_cfg.get('n_pcs', 30)
n_neighbors = pp_cfg.get('n_neighbors', 15)
leiden_resolution = pp_cfg.get('leiden_resolution', 0.5)
leiden_key = pp_cfg.get('leiden_key', 'leiden_res_0.50')

adata = sc.read_h5ad(src_h5ad)

# Basic QC metrics fallback-safe
if 'mt' not in adata.var.columns:
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
if 'ribo' not in adata.var.columns:
    adata.var['ribo'] = adata.var_names.str.upper().str.startswith(('RPS', 'RPL'))
if 'hb' not in adata.var.columns:
    adata.var['hb'] = adata.var_names.str.upper().str.startswith(('HBA', 'HBB'))

try:
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], percent_top=None, inplace=True)
except Exception:
    pass

# Optional filtering with config-defined cutoffs.
try:
    if min_genes_per_cell is not None and 'n_genes_by_counts' in adata.obs.columns:
        adata = adata[adata.obs['n_genes_by_counts'] >= float(min_genes_per_cell)].copy()
    if max_genes_per_cell is not None and 'n_genes_by_counts' in adata.obs.columns:
        adata = adata[adata.obs['n_genes_by_counts'] <= float(max_genes_per_cell)].copy()
    if min_counts_per_cell is not None and 'total_counts' in adata.obs.columns:
        adata = adata[adata.obs['total_counts'] >= float(min_counts_per_cell)].copy()
    if max_pct_counts_mt is not None and 'pct_counts_mt' in adata.obs.columns:
        adata = adata[adata.obs['pct_counts_mt'] <= float(max_pct_counts_mt)].copy()
    if min_cells_per_gene is not None and adata.n_vars > 0:
        sc.pp.filter_genes(adata, min_cells=max(1, int(min_cells_per_gene)))
except Exception:
    pass

# Compute PCA/neighbors/UMAP/Leiden on existing data if possible
try:
    if adata.n_vars > 1 and adata.n_obs > 2:
        if adata.X is not None:
            # use sparse-safe operations by scanpy
            sc.pp.normalize_total(adata, target_sum=float(target_sum))
            sc.pp.log1p(adata)
            n_hvg = min(int(n_top_genes), adata.n_vars)
            sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor=str(hvg_flavor), subset=False)
            sc.tl.pca(adata, svd_solver='arpack')
            n_neighbors_eff = min(int(n_neighbors), max(2, adata.n_obs - 1))
            n_pcs_eff = min(int(n_pcs), max(2, adata.n_vars - 1))
            sc.pp.neighbors(adata, n_neighbors=n_neighbors_eff, n_pcs=n_pcs_eff)
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=float(leiden_resolution), key_added=str(leiden_key))
except Exception:
    pass

# QC violin figure
qc_cols = [c for c in ['total_counts', 'n_genes_by_counts', 'pct_counts_mt'] if c in adata.obs.columns]
if qc_cols:
    sc.pl.violin(adata, qc_cols, jitter=0.4, multi_panel=True, show=False)
    ensure_parent(qc_png)
    plt.savefig(qc_png, dpi=200, bbox_inches='tight')
    plt.close('all')
else:
    save_placeholder(qc_png, 'QC metrics unavailable')

# UMAP figure
if 'X_umap' in adata.obsm:
    color_keys = [c for c in ['sample_id', 'final_label_transfer', str(leiden_key)] if c in adata.obs.columns]
    if not color_keys:
        color_keys = None
    sc.pl.umap(adata, color=color_keys, ncols=2, show=False)
    ensure_parent(umap_png)
    plt.savefig(umap_png, dpi=200, bbox_inches='tight')
    plt.close('all')
else:
    save_placeholder(umap_png, 'UMAP unavailable')

# Leiden figure
if 'X_umap' in adata.obsm and str(leiden_key) in adata.obs.columns:
    sc.pl.umap(adata, color=[str(leiden_key)], show=False)
    ensure_parent(leiden_png)
    plt.savefig(leiden_png, dpi=200, bbox_inches='tight')
    plt.close('all')
else:
    save_placeholder(leiden_png, 'Leiden unavailable')

ensure_parent(out_h5ad)
adata.write_h5ad(out_h5ad)
