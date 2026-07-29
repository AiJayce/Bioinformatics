import os
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc


def ensure_parent(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)


def save_placeholder(path, title):
    ensure_parent(path)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.text(0.5, 0.5, title, ha='center', va='center')
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


adata = sc.read_h5ad(snakemake.input.adata)
with open(snakemake.input.config, 'r', encoding='utf-8') as f:
    cfg = yaml.safe_load(f)

cfg_from_smk = getattr(snakemake, 'config', None)
if isinstance(cfg_from_smk, dict) and len(cfg_from_smk) > 0:
    cfg = cfg_from_smk

n_pcs_elbow = int(cfg.get('assignment1_notebook_figures', {}).get('n_pcs_elbow', 50))

out_pca = snakemake.output.pca_elbow
out_doublet = snakemake.output.doublet_hist
out_scvi = snakemake.output.scvi_loss
out_dotplot = snakemake.output.rank_dotplot
out_ccf_qc = snakemake.output.ccf_qc
out_bbknn_harmony = snakemake.output.bbknn_harmony
out_class_violin = snakemake.output.class_score_violin

# 1) PCA elbow
try:
    ad = adata.copy()
    if 'X_pca' not in ad.obsm or 'pca' not in ad.uns:
        sc.tl.pca(ad, n_comps=min(n_pcs_elbow, max(2, ad.n_vars - 1)))
    vr = np.asarray(ad.uns['pca']['variance_ratio'])
    ensure_parent(out_pca)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(np.arange(1, len(vr) + 1), vr, marker='o', linestyle='-')
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Explained Variance Ratio')
    ax.set_title('Elbow Plot for PCA')
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(out_pca, dpi=200)
    plt.close(fig)
except Exception as e:
    save_placeholder(out_pca, f'PCA elbow unavailable: {str(e)[:100]}')

# 2) Doublet score distribution
try:
    score_col = None
    for c in ['doublet_score', 'predicted_doublet', 'scrublet_score', 'scrublet_score_z']:
        if c in adata.obs.columns:
            score_col = c
            break
    ensure_parent(out_doublet)
    fig, ax = plt.subplots(figsize=(6, 4))
    if score_col is not None:
        vals = pd.to_numeric(adata.obs[score_col], errors='coerce').dropna().values
        ax.hist(vals, bins=50, alpha=0.85, color='#4C78A8')
        ax.set_title('Distribution of Doublet Scores')
        ax.set_xlabel(score_col)
    else:
        # fallback proxy if scrublet score not stored in adata
        vals = pd.to_numeric(adata.obs.get('n_genes_by_counts', pd.Series(dtype=float)), errors='coerce').dropna().values
        ax.hist(vals, bins=50, alpha=0.85, color='#4C78A8')
        ax.set_title('Distribution proxy (n_genes_by_counts)')
        ax.set_xlabel('n_genes_by_counts')
    ax.set_ylabel('Frequency')
    fig.tight_layout()
    fig.savefig(out_doublet, dpi=200)
    plt.close(fig)
except Exception as e:
    save_placeholder(out_doublet, f'Doublet histogram unavailable: {str(e)[:100]}')

# 3) scVI loss curve (from saved history if available, else placeholder)
try:
    ensure_parent(out_scvi)
    hist = adata.uns.get('scvi_history', None)
    fig, ax = plt.subplots(figsize=(7, 4))
    if isinstance(hist, dict) and 'elbo_train' in hist:
        train = np.asarray(hist.get('elbo_train', []), dtype=float)
        val = np.asarray(hist.get('elbo_validation', []), dtype=float)
        if len(train) > 0:
            ax.plot(train, label='train')
        if len(val) > 0:
            ax.plot(val, label='validation')
        ax.legend()
        ax.set_title('SCVI Training and Validation Loss')
    else:
        ax.text(0.5, 0.5, 'scVI history unavailable in adata.uns', ha='center', va='center')
        ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_scvi, dpi=200)
    plt.close(fig)
except Exception as e:
    save_placeholder(out_scvi, f'scVI loss unavailable: {str(e)[:100]}')

# 4) rank_genes_groups dotplot
try:
    if 'leiden_res_0.25' in adata.obs.columns:
        group_col = 'leiden_res_0.25'
    elif 'leiden_res_0.10' in adata.obs.columns:
        group_col = 'leiden_res_0.10'
    elif 'cluster_pred' in adata.obs.columns:
        group_col = 'cluster_pred'
    elif 'final_label_transfer' in adata.obs.columns:
        group_col = 'final_label_transfer'
    elif 'mannual_annotaion' in adata.obs.columns:
        group_col = 'mannual_annotaion'
    else:
        group_col = None
    if group_col is None:
        raise ValueError('No group column for rank_genes_groups')

    ad = adata.copy()
    ad.obs[group_col] = ad.obs[group_col].astype(str).astype('category')
    sc.tl.rank_genes_groups(ad, groupby=group_col, method='wilcoxon', use_raw=False)
    sc.pl.rank_genes_groups_dotplot(ad, groupby=group_col, n_genes=5, show=False)
    ensure_parent(out_dotplot)
    plt.savefig(out_dotplot, dpi=200, bbox_inches='tight')
    plt.close('all')
except Exception as e:
    save_placeholder(out_dotplot, f'rank_genes_groups dotplot unavailable: {str(e)[:100]}')

# 5) CCF registration/QC panel
try:
    req = ['center_x', 'center_y', 'x_mosaic_px', 'y_mosaic_px', 'CCF_acronym']
    if any(c not in adata.obs.columns for c in req):
        raise ValueError('Missing CCF/QC columns')
    ensure_parent(out_ccf_qc)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    axes[0].scatter(adata.obs['center_x'], adata.obs['center_y'], s=1, alpha=0.15)
    axes[0].set_title('Original coordinates')
    axes[0].invert_yaxis()

    axes[1].scatter(adata.obs['x_mosaic_px'], adata.obs['y_mosaic_px'], s=1, alpha=0.15)
    axes[1].set_title('Mosaic pixel coordinates')
    axes[1].invert_yaxis()

    top = adata.obs['CCF_acronym'].astype(str).value_counts().head(20)
    top.plot(kind='bar', ax=axes[2])
    axes[2].set_title('Top CCF labels')
    axes[2].tick_params(axis='x', rotation=70)

    fig.tight_layout()
    fig.savefig(out_ccf_qc, dpi=200)
    plt.close(fig)
except Exception as e:
    save_placeholder(out_ccf_qc, f'CCF QC panel unavailable: {str(e)[:100]}')

# 6) BBKNN/Harmony UMAP compare (best-effort)
try:
    ensure_parent(out_bbknn_harmony)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    if 'X_umap' in adata.obsm:
        xy = np.asarray(adata.obsm['X_umap'])
        axes[0].scatter(xy[:, 0], xy[:, 1], s=2, alpha=0.3)
        axes[0].set_title('UMAP (reference)')
        axes[1].scatter(xy[:, 0], xy[:, 1], s=2, alpha=0.3)
        axes[1].set_title('UMAP (BBKNN/Harmony placeholder)')
    else:
        for ax in axes:
            ax.text(0.5, 0.5, 'UMAP unavailable', ha='center', va='center')
            ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_bbknn_harmony, dpi=200)
    plt.close(fig)
except Exception as e:
    save_placeholder(out_bbknn_harmony, f'BBKNN/Harmony compare unavailable: {str(e)[:100]}')

# 7) class_score violin (best-effort)
try:
    score_cols = [c for c in adata.obs.columns if c.endswith('_score') or c in ['class_score', 'CCF_score']]
    if len(score_cols) == 0:
        raise ValueError('No score columns found')
    cols = score_cols[:6]
    sc.pl.violin(adata, cols, jitter=0.2, multi_panel=True, show=False)
    ensure_parent(out_class_violin)
    plt.savefig(out_class_violin, dpi=200, bbox_inches='tight')
    plt.close('all')
except Exception as e:
    save_placeholder(out_class_violin, f'Class-score violin unavailable: {str(e)[:100]}')
