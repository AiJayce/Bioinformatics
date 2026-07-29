import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import scanpy as sc


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


latent_csv = snakemake.input.latent
adata_h5ad = snakemake.input.adata
micro_vsd_csv = snakemake.input.micro_vsd
oligo_vsd_csv = snakemake.input.oligo_vsd
micro_pseudobulk_csv = snakemake.input.micro_pseudobulk
oligo_pseudobulk_csv = snakemake.input.oligo_pseudobulk

latent_dim_count_png = snakemake.output.latent_dim_count
effect_dist_png = snakemake.output.effect_dist
latent_heatmap_png = snakemake.output.latent_heatmap
z15_bar_png = snakemake.output.z15_bar
micro_pca_png = snakemake.output.micro_pca
oligo_pca_png = snakemake.output.oligo_pca
neighbor_heatmap_png = snakemake.output.neighbor_heatmap
trem2_scatter_png = snakemake.output.trem2_scatter
apoe_scatter_png = snakemake.output.apoe_scatter
micro_pseudobulk_heatmap_png = snakemake.output.micro_pseudobulk_heatmap
oligo_pseudobulk_heatmap_png = snakemake.output.oligo_pseudobulk_heatmap

# 1) Latent-based plots
df = pd.read_csv(latent_csv)
df.columns = [c.strip().lstrip('\ufeff') for c in df.columns]
required = {'gene', 'latent_dim', 'effect'}
if required.issubset(df.columns):
    df = df.copy()
    df['gene'] = df['gene'].astype(str)
    df['latent_dim'] = df['latent_dim'].astype(str)
    df['effect'] = pd.to_numeric(df['effect'], errors='coerce')
    df = df.dropna(subset=['effect'])

    # Plot A: genes per latent dimension
    ensure_parent(latent_dim_count_png)
    cnt = df.groupby('latent_dim')['gene'].nunique().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(10, 4))
    cnt.plot(kind='bar', ax=ax, color='#4C78A8')
    ax.set_title('Assignment_2: Gene Count per Latent Dimension')
    ax.set_xlabel('latent_dim')
    ax.set_ylabel('unique gene count')
    plt.tight_layout()
    fig.savefig(latent_dim_count_png, dpi=200)
    plt.close(fig)

    # Plot B: effect distribution
    ensure_parent(effect_dist_png)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(df['effect'].values, bins=60, color='#F58518', alpha=0.85)
    ax.set_title('Assignment_2: Latent Effect Distribution')
    ax.set_xlabel('effect')
    ax.set_ylabel('frequency')
    plt.tight_layout()
    fig.savefig(effect_dist_png, dpi=200)
    plt.close(fig)

    # Plot C: heatmap (top abs effect genes per latent)
    ensure_parent(latent_heatmap_png)
    top_n = 15
    top_rows = (
        df.assign(abs_effect=df['effect'].abs())
          .sort_values(['latent_dim', 'abs_effect'], ascending=[True, False])
          .groupby('latent_dim')
          .head(top_n)
    )
    heat = top_rows.pivot_table(index='gene', columns='latent_dim', values='effect', aggfunc='mean').fillna(0)
    if heat.shape[0] > 0 and heat.shape[1] > 0:
        fig, ax = plt.subplots(figsize=(max(8, 0.45 * heat.shape[1]), max(6, 0.18 * heat.shape[0])))
        im = ax.imshow(heat.values, aspect='auto', cmap='coolwarm', interpolation='nearest')
        ax.set_xticks(np.arange(heat.shape[1]))
        ax.set_xticklabels(heat.columns, rotation=90)
        ax.set_yticks(np.arange(heat.shape[0]))
        ax.set_yticklabels(heat.index)
        ax.set_title('Assignment_2: Top Latent Genes Heatmap')
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label('effect')
        plt.tight_layout()
        fig.savefig(latent_heatmap_png, dpi=200)
        plt.close(fig)
    else:
        save_placeholder(latent_heatmap_png, 'No data for latent heatmap')

    # Plot D: z15 top genes
    ensure_parent(z15_bar_png)
    z15 = df[df['latent_dim'] == 'z15'].copy()
    if len(z15) == 0:
        z15 = df.copy()
    z15 = z15.assign(abs_effect=z15['effect'].abs()).sort_values('abs_effect', ascending=False).head(20)
    if len(z15) > 0:
        z15 = z15.sort_values('effect')
        fig, ax = plt.subplots(figsize=(8, 6))
        colors = np.where(z15['effect'].values >= 0, '#54A24B', '#E45756')
        ax.barh(z15['gene'].values, z15['effect'].values, color=colors)
        ax.set_title('Assignment_2: Top Genes in z15 (or global fallback)')
        ax.set_xlabel('effect')
        ax.set_ylabel('gene')
        plt.tight_layout()
        fig.savefig(z15_bar_png, dpi=200)
        plt.close(fig)
    else:
        save_placeholder(z15_bar_png, 'No data for z15 top genes')
else:
    save_placeholder(latent_dim_count_png, 'Missing columns in latent CSV')
    save_placeholder(effect_dist_png, 'Missing columns in latent CSV')
    save_placeholder(latent_heatmap_png, 'Missing columns in latent CSV')
    save_placeholder(z15_bar_png, 'Missing columns in latent CSV')


# 2) DESeq2 VSD PCA plots
for in_csv, out_png, title in [
    (micro_vsd_csv, micro_pca_png, 'Assignment_2: Microglia VSD PCA'),
    (oligo_vsd_csv, oligo_pca_png, 'Assignment_2: Oligodendrocyte VSD PCA'),
]:
    expr = pd.read_csv(in_csv)
    expr.columns = [c.strip().lstrip('\ufeff') for c in expr.columns]
    if expr.shape[1] < 3:
        save_placeholder(out_png, f'Not enough columns for PCA: {os.path.basename(in_csv)}')
        continue

    first_col = expr.columns[0]
    mat = expr.set_index(first_col)
    mat = mat.apply(pd.to_numeric, errors='coerce').fillna(0.0)

    # samples x genes
    X = mat.T.values
    if X.shape[0] < 2 or X.shape[1] < 2:
        save_placeholder(out_png, f'Not enough matrix shape for PCA: {os.path.basename(in_csv)}')
        continue

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)

    ensure_parent(out_png)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(pcs[:, 0], pcs[:, 1], s=30, alpha=0.85, color='#4C78A8')
    for i, sample in enumerate(mat.columns):
        ax.text(pcs[i, 0], pcs[i, 1], str(sample), fontsize=7, alpha=0.8)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)')
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def pick_label_column(obs):
    for c in ['cluster_pred', 'final_label_transfer', 'mannual_annotaion', 'Mannual_annotaion', 'class']:
        if c in obs.columns:
            return c
    return None


def to_dense_vector(x):
    if hasattr(x, 'toarray'):
        x = x.toarray()
    return np.asarray(x).reshape(-1)


# 3) Microglia Trem2/Apoe Neighbor Association (MERSCOPE-like)
try:
    adata = sc.read_h5ad(adata_h5ad)
    obs = adata.obs.copy()
    label_col = pick_label_column(obs)

    if label_col is None or 'sample_id' not in obs.columns:
        raise ValueError('Required columns for neighbor analysis are missing')

    if 'center_x' in obs.columns and 'center_y' in obs.columns:
        coord_source = 'obs'
    elif hasattr(adata, 'obsm') and 'spatial' in adata.obsm:
        coord_source = 'obsm_spatial'
    else:
        raise ValueError('No usable spatial coordinates (center_x/center_y or obsm[spatial])')

    label_series = obs[label_col].astype(str)
    micro_mask = label_series.str.contains('Microglia', case=False, na=False)
    if micro_mask.sum() < 3:
        raise ValueError('Not enough Microglia cells for neighbor analysis')

    genes_upper = pd.Index([g.upper() for g in adata.var_names])
    if 'TREM2' not in genes_upper or 'APOE' not in genes_upper:
        raise ValueError('Trem2/Apoe genes not found in adata.var_names')

    trem2_idx = int(np.where(genes_upper == 'TREM2')[0][0])
    apoe_idx = int(np.where(genes_upper == 'APOE')[0][0])

    micro_df_rows = []
    n_neighbors = 12

    for sid in obs['sample_id'].astype(str).unique():
        smask = obs['sample_id'].astype(str) == sid
        sub_obs = obs.loc[smask].copy()
        sub_idx = np.where(smask.values)[0]
        if len(sub_obs) < 5:
            continue

        if coord_source == 'obs':
            coords = sub_obs[['center_x', 'center_y']].to_numpy(dtype=float)
        else:
            coords = np.asarray(adata.obsm['spatial'])[sub_idx, :2]
        nbr_n = min(n_neighbors + 1, len(sub_obs))
        if nbr_n <= 2:
            continue

        nn = NearestNeighbors(n_neighbors=nbr_n)
        nn.fit(coords)
        _, idx = nn.kneighbors(coords)
        idx = idx[:, 1:]

        sub_labels = sub_obs[label_col].astype(str).values
        sub_micro = sub_labels.copy()
        sub_gene_trem2 = to_dense_vector(adata.X[sub_idx, trem2_idx])
        sub_gene_apoe = to_dense_vector(adata.X[sub_idx, apoe_idx])

        for i_local, global_i in enumerate(sub_idx):
            if 'microglia' not in sub_labels[i_local].lower():
                continue
            nbr = idx[i_local]
            nbr_labels = sub_labels[nbr]
            total = len(nbr_labels)
            if total == 0:
                continue

            uniq, counts = np.unique(nbr_labels, return_counts=True)
            row = {
                'sample_id': sid,
                'cell_index': int(global_i),
                'log1p_Trem2': float(np.log1p(sub_gene_trem2[i_local])),
                'log1p_Apoe': float(np.log1p(sub_gene_apoe[i_local])),
                'n_neighbors': int(total),
            }
            for u, c in zip(uniq, counts):
                key = f'nbr_frac::{u}'
                row[key] = float(c / total)
            micro_df_rows.append(row)

    micro_df = pd.DataFrame(micro_df_rows)
    frac_cols = sorted([c for c in micro_df.columns if c.startswith('nbr_frac::')])
    if len(micro_df) < 5 or len(frac_cols) == 0:
        raise ValueError('Insufficient neighbor table for association plots')

    # Keep top partner labels by average fraction.
    top_cols = (
        micro_df[frac_cols]
        .mean(axis=0)
        .sort_values(ascending=False)
        .head(12)
        .index
        .tolist()
    )

    assoc_rows = []
    for response_col, response_name in [('log1p_Trem2', 'Trem2'), ('log1p_Apoe', 'Apoe')]:
        for fc in top_cols:
            x = micro_df[fc].to_numpy(dtype=float)
            y = micro_df[response_col].to_numpy(dtype=float)
            if np.std(x) == 0 or np.std(y) == 0:
                corr = 0.0
            else:
                corr = float(np.corrcoef(x, y)[0, 1])
            assoc_rows.append({'gene': response_name, 'partner': fc.replace('nbr_frac::', ''), 'corr': corr})
    assoc_df = pd.DataFrame(assoc_rows)

    # Heatmap figure
    hm = assoc_df.pivot(index='gene', columns='partner', values='corr').fillna(0.0)
    ensure_parent(neighbor_heatmap_png)
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * hm.shape[1]), 4.5))
    im = ax.imshow(hm.values, aspect='auto', cmap='coolwarm', vmin=-1, vmax=1)
    ax.set_xticks(np.arange(hm.shape[1]))
    ax.set_xticklabels(hm.columns, rotation=60, ha='right')
    ax.set_yticks(np.arange(hm.shape[0]))
    ax.set_yticklabels(hm.index)
    ax.set_title('Microglia neighbor association with Trem2/Apoe')
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('Pearson r')
    plt.tight_layout()
    fig.savefig(neighbor_heatmap_png, dpi=200)
    plt.close(fig)

    # Scatter figures for Oligodendrocyte neighbors
    oligo_candidates = [c for c in frac_cols if 'oligodendro' in c.lower()]
    if len(oligo_candidates) == 0:
        oligo_col = top_cols[0]
    else:
        oligo_col = oligo_candidates[0]

    for ycol, out_png, title in [
        ('log1p_Trem2', trem2_scatter_png, 'Microglia Trem2 vs Oligodendrocyte neighbor fraction (k=12)'),
        ('log1p_Apoe', apoe_scatter_png, 'Microglia Apoe vs Oligodendrocyte neighbor fraction (k=12)'),
    ]:
        x = micro_df[oligo_col].to_numpy(dtype=float)
        y = micro_df[ycol].to_numpy(dtype=float)
        ensure_parent(out_png)
        fig, ax = plt.subplots(figsize=(6.8, 5.2))
        ax.scatter(x, y, s=8, alpha=0.35, color='#4C78A8')
        if len(x) > 3 and np.std(x) > 0:
            coef = np.polyfit(x, y, 1)
            xs = np.linspace(float(np.min(x)), float(np.max(x)), 100)
            ax.plot(xs, coef[0] * xs + coef[1], color='#E45756', linewidth=2)
        r = float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else 0.0
        ax.set_title(title)
        ax.set_xlabel(oligo_col.replace('nbr_frac::', '') + ' neighbor fraction')
        ax.set_ylabel(ycol)
        ax.text(0.02, 0.98, f'r = {r:.3f}', transform=ax.transAxes, va='top', ha='left')
        plt.tight_layout()
        fig.savefig(out_png, dpi=200)
        plt.close(fig)

except Exception as e:
    save_placeholder(neighbor_heatmap_png, f'Neighbor association unavailable: {str(e)[:80]}')
    save_placeholder(trem2_scatter_png, f'Trem2 scatter unavailable: {str(e)[:80]}')
    save_placeholder(apoe_scatter_png, f'Apoe scatter unavailable: {str(e)[:80]}')


def save_pseudobulk_heatmap(in_csv, out_png, title):
    expr = pd.read_csv(in_csv)
    expr.columns = [c.strip().lstrip('\ufeff') for c in expr.columns]
    if expr.shape[1] < 3:
        save_placeholder(out_png, f'Not enough columns: {os.path.basename(in_csv)}')
        return

    first_col = expr.columns[0]
    mat = expr.set_index(first_col)
    mat = mat.apply(pd.to_numeric, errors='coerce').fillna(0.0)
    if mat.shape[0] < 2 or mat.shape[1] < 2:
        save_placeholder(out_png, f'Not enough matrix shape: {os.path.basename(in_csv)}')
        return

    # top variable genes
    var = mat.var(axis=1).sort_values(ascending=False)
    top = var.head(min(50, len(var))).index
    sub = mat.loc[top]

    # z-score by gene
    mu = sub.mean(axis=1)
    sd = sub.std(axis=1).replace(0, 1.0)
    z = sub.sub(mu, axis=0).div(sd, axis=0)

    ensure_parent(out_png)
    fig, ax = plt.subplots(figsize=(max(7, 0.3 * z.shape[1]), max(6, 0.15 * z.shape[0])))
    im = ax.imshow(z.values, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
    ax.set_xticks(np.arange(z.shape[1]))
    ax.set_xticklabels(z.columns, rotation=90)
    ax.set_yticks(np.arange(z.shape[0]))
    ax.set_yticklabels(z.index)
    ax.set_title(title)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('gene-wise z-score')
    plt.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# 4) NicheNet pseudobulk-style heatmaps
save_pseudobulk_heatmap(micro_pseudobulk_csv, micro_pseudobulk_heatmap_png, 'NicheNet/Pseudobulk: Microglia top variable genes')
save_pseudobulk_heatmap(oligo_pseudobulk_csv, oligo_pseudobulk_heatmap_png, 'NicheNet/Pseudobulk: Oligodendrocyte top variable genes')
