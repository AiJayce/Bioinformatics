import os
import glob
import re
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import sparse


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


def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def to_dense(x):
    if sparse.issparse(x):
        return x.toarray()
    if hasattr(x, 'toarray'):
        return x.toarray()
    return np.asarray(x)


adata_path = snakemake.input.adata
latent_path = snakemake.input.latent
vae_model_dir = snakemake.input.vae_model_dir
plaques_h5ad_path = snakemake.input.plaques_h5ad
oligo_vsd_path = snakemake.input.oligo_vsd
micro_vsd_path = snakemake.input.micro_vsd

out_boxplot = snakemake.output.boxplot
out_volcano_wt = snakemake.output.volcano_wt
out_volcano_ad = snakemake.output.volcano_ad
out_venn = snakemake.output.venn
out_spatial = snakemake.output.spatial
out_spatial_zoom = snakemake.output.spatial_zoom
out_vae_heatmap = snakemake.output.vae_heatmap
out_vae_top10 = snakemake.output.vae_top10
out_vae_dist = snakemake.output.vae_distribution
out_plaque_overlay = snakemake.output.plaque_overlay
out_plaque_bar = snakemake.output.plaque_bar
out_lr_scatter = snakemake.output.lr_scatter
out_umap_stab1_ntm = snakemake.output.umap_stab1_ntm
out_umap_manual = snakemake.output.umap_manual

with open(snakemake.input.config, 'r', encoding='utf-8') as f:
    cfg = yaml.safe_load(f)

cfg_from_smk = getattr(snakemake, 'config', None)
if isinstance(cfg_from_smk, dict) and len(cfg_from_smk) > 0:
    cfg = cfg_from_smk

fig_cfg = cfg.get('assignment2_notebook_figures', {})
micro_label = fig_cfg.get('microglia_label', 'Microglia')
oligo_label = fig_cfg.get('oligodendrocyte_label', 'Oligodendrocyte')
posterior_label = fig_cfg.get('posterior_label', 'Posterior Hippocampus')
region_major_label = fig_cfg.get('region_major_label', 'Hippocampus')
vae_focus_latent = str(fig_cfg.get('vae_focus_latent', 'z15'))
vae_top_n_genes = int(fig_cfg.get('vae_top_n_genes', 10))
plaque_sample_map_cfg = fig_cfg.get('plaque_sample_map', [])

adata = sc.read_h5ad(adata_path)
obs = adata.obs.copy()
var_upper = pd.Index([g.upper() for g in adata.var_names])

sample_col = pick_col(obs, ['sample_id', 'library_label'])
celltype_col = pick_col(obs, ['cluster_pred', 'final_label_transfer', 'mannual_annotaion', 'Mannual_annotaion', 'class'])
genotype_col = pick_col(obs, ['genotype', 'donor_genotype'])
region_col = pick_col(obs, ['region_location', 'region_of_interest_acronym'])
major_col = pick_col(obs, ['CCF_major'])

# 1) Trem2/Apoe boxplot by genotype (and optionally region)
try:
    if celltype_col is None or genotype_col is None:
        raise ValueError('missing cell type or genotype column')
    micro = adata[obs[celltype_col].astype(str).str.contains(micro_label, case=False, na=False)].copy()
    if micro.n_obs < 5:
        raise ValueError('too few microglia cells')

    genes = []
    for g in ['TREM2', 'APOE']:
        if g in var_upper:
            genes.append(adata.var_names[int(np.where(var_upper == g)[0][0])])
    if len(genes) < 2:
        raise ValueError('Trem2/Apoe not found')

    X = to_dense(micro[:, genes].X)
    plot_df = pd.DataFrame(X, columns=genes, index=micro.obs_names)
    plot_df['genotype'] = micro.obs[genotype_col].astype(str).replace({'5xFAD': 'AD'})
    if region_col is not None:
        plot_df['region'] = micro.obs[region_col].astype(str)
    m = plot_df.melt(id_vars=[c for c in ['genotype', 'region'] if c in plot_df.columns], var_name='gene', value_name='expression')

    ensure_parent(out_boxplot)
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.boxplot(data=m, x='gene', y='expression', hue='genotype', ax=ax, showfliers=False)
    ax.set_title('Microglia Trem2/Apoe expression by genotype')
    ax.set_xlabel('')
    ax.set_ylabel('Expression')
    plt.tight_layout()
    fig.savefig(out_boxplot, dpi=200)
    plt.close(fig)
except Exception as e:
    save_placeholder(out_boxplot, f'Boxplot unavailable: {str(e)[:100]}')


def make_volcano(sub, out_path, title):
    try:
        if region_col is None:
            raise ValueError('region column missing')
        if sub.n_obs < 20:
            raise ValueError('too few cells for DEG')
        groups = sub.obs[region_col].astype(str)
        if posterior_label not in set(groups):
            raise ValueError('posterior label missing')
        sc.tl.rank_genes_groups(sub, groupby=region_col, groups=[posterior_label], reference='rest', method='wilcoxon')
        deg = sc.get.rank_genes_groups_df(sub, group=posterior_label)
        deg['neglog10'] = -np.log10(np.clip(deg['pvals'].astype(float), 1e-300, 1.0))
        deg['signif'] = 'NS'
        deg.loc[(deg['pvals'] < 0.05) & (deg['logfoldchanges'] > 0.5), 'signif'] = 'Up'
        deg.loc[(deg['pvals'] < 0.05) & (deg['logfoldchanges'] < -0.5), 'signif'] = 'Down'

        ensure_parent(out_path)
        fig, ax = plt.subplots(figsize=(6, 5))
        cmap = {'NS': '#999999', 'Up': '#D62728', 'Down': '#1F77B4'}
        for k in ['NS', 'Down', 'Up']:
            d = deg[deg['signif'] == k]
            ax.scatter(d['logfoldchanges'], d['neglog10'], s=8, c=cmap[k], alpha=0.7)
        ax.axvline(0, color='black', linewidth=1)
        ax.axvline(0.5, color='black', linewidth=1, linestyle='--')
        ax.axvline(-0.5, color='black', linewidth=1, linestyle='--')
        ax.axhline(-np.log10(0.05), color='black', linewidth=1, linestyle='--')
        ax.set_xlabel('log2 Fold Change (Posterior vs rest)')
        ax.set_ylabel('-log10(p-value)')
        ax.set_title(title)
        plt.tight_layout()
        fig.savefig(out_path, dpi=200)
        plt.close(fig)
        return set(deg.loc[(deg['pvals'] < 0.05) & (deg['logfoldchanges'] > 0.5), 'names'].astype(str))
    except Exception as e:
        save_placeholder(out_path, f'Volcano unavailable: {str(e)[:100]}')
        return set()


# 2) DEG volcano WT/AD and Venn
wt_up = set()
ad_up = set()
try:
    if celltype_col is None or genotype_col is None:
        raise ValueError('missing required columns')
    sub = adata.copy()
    if major_col is not None:
        sub = sub[sub.obs[major_col].astype(str) == region_major_label].copy()
    sub = sub[sub.obs[celltype_col].astype(str).str.contains(micro_label, case=False, na=False)].copy()

    wt = sub[sub.obs[genotype_col].astype(str).isin(['WT'])].copy()
    ad = sub[sub.obs[genotype_col].astype(str).isin(['5xFAD', 'AD'])].copy()

    wt_up = make_volcano(wt, out_volcano_wt, 'WT: Posterior vs Anterior/Rest')
    ad_up = make_volcano(ad, out_volcano_ad, 'AD: Posterior vs Anterior/Rest')
except Exception as e:
    save_placeholder(out_volcano_wt, f'WT volcano unavailable: {str(e)[:100]}')
    save_placeholder(out_volcano_ad, f'AD volcano unavailable: {str(e)[:100]}')

try:
    ensure_parent(out_venn)
    fig, ax = plt.subplots(figsize=(7, 6))
    try:
        from matplotlib_venn import venn2
        v = venn2([wt_up, ad_up], set_labels=('WT up', 'AD up'), ax=ax)
        if v is None:
            raise ValueError('venn2 failed')
        ax.set_title('Posterior-up DEG overlap (Microglia)')
    except Exception:
        inter = len(wt_up & ad_up)
        ax.text(0.5, 0.65, f'WT up: {len(wt_up)}', ha='center')
        ax.text(0.5, 0.50, f'AD up: {len(ad_up)}', ha='center')
        ax.text(0.5, 0.35, f'Shared: {inter}', ha='center', fontweight='bold')
        ax.set_axis_off()
    plt.tight_layout()
    fig.savefig(out_venn, dpi=200)
    plt.close(fig)
except Exception as e:
    save_placeholder(out_venn, f'VENN unavailable: {str(e)[:100]}')


# 3) Spatial overlay Apoe/Trem2 in Microglia + Oligo outline
try:
    if celltype_col is None:
        raise ValueError('missing cell type column')
    sub = adata[adata.obs[celltype_col].astype(str).isin([micro_label, oligo_label])].copy()
    if sub.n_obs < 10:
        raise ValueError('too few cells for spatial overlay')

    if 'spatial' in sub.obsm:
        xy = np.asarray(sub.obsm['spatial'])
    else:
        if 'center_x' not in sub.obs.columns or 'center_y' not in sub.obs.columns:
            raise ValueError('missing spatial coordinates')
        xy = sub.obs[['center_x', 'center_y']].to_numpy(dtype=float)

    g1 = adata.var_names[int(np.where(var_upper == 'APOE')[0][0])] if 'APOE' in var_upper else None
    g2 = adata.var_names[int(np.where(var_upper == 'TREM2')[0][0])] if 'TREM2' in var_upper else None
    if g1 is None or g2 is None:
        raise ValueError('Apoe/Trem2 genes missing')

    micro_mask = sub.obs[celltype_col].astype(str).str.lower() == micro_label.lower()
    oligo_mask = sub.obs[celltype_col].astype(str).str.lower() == oligo_label.lower()
    if micro_mask.sum() < 2 or oligo_mask.sum() < 2:
        raise ValueError('insufficient microglia/oligodendrocyte cells')

    Xg = to_dense(sub[:, [g1, g2]].X)
    expr = np.mean(Xg, axis=1)
    expr_min = float(np.min(expr))
    expr_max = float(np.max(expr))
    if expr_max > expr_min:
        expr = (expr - expr_min) / (expr_max - expr_min)
    else:
        expr = np.zeros_like(expr)

    def draw_spatial(out_path, zoom=False):
        ensure_parent(out_path)
        fig, ax = plt.subplots(figsize=(7.5, 6.5))
        ax.scatter(xy[oligo_mask, 0], xy[oligo_mask, 1], s=8 if not zoom else 18, facecolors='none', edgecolors='blue', linewidths=0.5)
        sca = ax.scatter(xy[micro_mask, 0], xy[micro_mask, 1], s=8 if not zoom else 22, c=expr[micro_mask], cmap='Greens', edgecolors='red', linewidths=0.5)
        if zoom:
            mx = xy[micro_mask, 0]
            my = xy[micro_mask, 1]
            cx, cy = np.mean(mx), np.mean(my)
            r = max(np.std(mx), np.std(my)) * 2.0
            if r <= 0:
                r = 500
            ax.set_xlim(cx - r, cx + r)
            ax.set_ylim(cy + r, cy - r)
        else:
            ax.invert_yaxis()
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('Microglia Apoe/Trem2 with Oligodendrocyte context' + (' (zoom)' if zoom else ''))
        fig.colorbar(sca, ax=ax, fraction=0.05, pad=0.04, label='mean scaled expr (Apoe/Trem2)')
        plt.tight_layout()
        fig.savefig(out_path, dpi=200)
        plt.close(fig)

    draw_spatial(out_spatial, zoom=False)
    draw_spatial(out_spatial_zoom, zoom=True)
except Exception as e:
    save_placeholder(out_spatial, f'Spatial overlay unavailable: {str(e)[:100]}')
    save_placeholder(out_spatial_zoom, f'Spatial zoom unavailable: {str(e)[:100]}')


# 4) VAE reproducibility figures from saved ridge coefficients and latent effects
latent_df = pd.read_csv(latent_path)
latent_df.columns = [c.strip().lstrip('\ufeff') for c in latent_df.columns]
if 'latent_dim' in latent_df.columns:
    latent_df['latent_dim'] = latent_df['latent_dim'].astype(str)
if 'effect' in latent_df.columns:
    latent_df['effect'] = pd.to_numeric(latent_df['effect'], errors='coerce')

# 4a) Heatmap from all celltype ridge_coefficient.csv if available
try:
    ridge_files = sorted(glob.glob(os.path.join(vae_model_dir, '*', 'ridge_coefficient.csv')))
    mats = []
    for fp in ridge_files:
        ct = os.path.basename(os.path.dirname(fp))
        d = pd.read_csv(fp)
        d.columns = [c.strip().lstrip('\ufeff') for c in d.columns]
        if 'latent_dim' in d.columns and 'coefficient' in d.columns:
            p = d.pivot_table(index=None, columns='latent_dim', values='coefficient', aggfunc='mean')
            p.index = [ct]
            mats.append(p)
    if len(mats) == 0:
        raise ValueError('no ridge coefficient files found')
    hm = pd.concat(mats, axis=0).fillna(0.0)
    hm = hm.reindex(sorted(hm.columns), axis=1)
    ensure_parent(out_vae_heatmap)
    fig, ax = plt.subplots(figsize=(max(8, 0.6 * hm.shape[1]), max(4, 0.5 * hm.shape[0])))
    im = ax.imshow(hm.values, aspect='auto', cmap='RdBu_r', vmin=-0.5, vmax=0.5)
    ax.set_xticks(np.arange(hm.shape[1]))
    ax.set_xticklabels(hm.columns, rotation=60, ha='right')
    ax.set_yticks(np.arange(hm.shape[0]))
    ax.set_yticklabels(hm.index)
    ax.set_title('VAE latent coefficients across cell types')
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('coefficient')
    plt.tight_layout()
    fig.savefig(out_vae_heatmap, dpi=200)
    plt.close(fig)
except Exception as e:
    save_placeholder(out_vae_heatmap, f'VAE heatmap unavailable: {str(e)[:100]}')

# 4b) z15 (or fallback) top gene contribution bar
try:
    req = {'gene', 'latent_dim', 'effect'}
    if not req.issubset(latent_df.columns):
        raise ValueError('latent file missing required columns')
    z = latent_df[latent_df['latent_dim'] == vae_focus_latent].copy()
    if len(z) == 0:
        z = latent_df.copy()
    z = z.dropna(subset=['effect']).sort_values('effect', ascending=False).head(vae_top_n_genes).sort_values('effect')
    if len(z) == 0:
        raise ValueError('no rows for top gene plot')
    ensure_parent(out_vae_top10)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.barh(z['gene'].astype(str), z['effect'].astype(float), color=np.where(z['effect'].astype(float) >= 0, '#54A24B', '#E45756'))
    ax.axvline(0, color='black', linewidth=1)
    ax.set_title(f'Top {vae_top_n_genes} genes in {vae_focus_latent} (fallback enabled)')
    ax.set_xlabel('effect')
    ax.set_ylabel('')
    plt.tight_layout()
    fig.savefig(out_vae_top10, dpi=200)
    plt.close(fig)
except Exception as e:
    save_placeholder(out_vae_top10, f'VAE top genes unavailable: {str(e)[:100]}')

# 4c) Distribution of top genes in microglia
try:
    req = {'gene', 'latent_dim', 'effect'}
    if not req.issubset(latent_df.columns):
        raise ValueError('latent file missing required columns')
    top_genes = (
        latent_df.assign(abs_effect=latent_df['effect'].abs())
        .sort_values('abs_effect', ascending=False)
        .head(vae_top_n_genes)['gene']
        .astype(str)
        .tolist()
    )
    top_genes = [g for g in top_genes if g in adata.var_names]
    if len(top_genes) == 0:
        raise ValueError('top genes not found in adata.var_names')
    if celltype_col is None:
        raise ValueError('missing cell type column')
    micro = adata[adata.obs[celltype_col].astype(str).str.contains(micro_label, case=False, na=False)].copy()
    if micro.n_obs < 5:
        raise ValueError('too few microglia cells')

    X = to_dense(micro[:, top_genes].X)
    expr = pd.DataFrame(X, columns=top_genes)
    ensure_parent(out_vae_dist)
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    for g in top_genes:
        sns.histplot(expr[g], bins=80, stat='density', element='step', fill=False, linewidth=2, ax=ax, label=g)
    ax.set_title(f'Microglia expression distribution of top {len(top_genes)} VAE genes')
    ax.set_xlabel('Expression')
    ax.set_ylabel('Density')
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), frameon=True)
    plt.tight_layout()
    fig.savefig(out_vae_dist, dpi=200)
    plt.close(fig)
except Exception as e:
    save_placeholder(out_vae_dist, f'VAE distribution unavailable: {str(e)[:100]}')


# 5) Plaque + MERFISH overlay and plaque burden bar
try:
    pdata = sc.read_h5ad(plaques_h5ad_path)

    def make_key(path):
        p = str(path).replace('\\\\', '/').replace('\\', '/')
        m = re.search(r'region_(\d+)', p)
        if m is None:
            return None
        region = m.group(1)
        sample_part = p.split('/region_')[0].rstrip('/')
        sample = sample_part.split('/')[-1] if sample_part else ''
        if not sample:
            return None
        return f'{sample}_region_{region}'

    if 'sample' in pdata.obs.columns:
        pdata.obs['key'] = pdata.obs['sample'].astype(str).apply(make_key)
    else:
        pdata.obs['key'] = None

    if not isinstance(plaque_sample_map_cfg, list) or len(plaque_sample_map_cfg) == 0:
        raise ValueError('assignment2_notebook_figures.plaque_sample_map is missing')

    data_info = pd.DataFrame(plaque_sample_map_cfg)
    required_cols = {'sample_path', 'sample_id', 'genotype', 'region_location'}
    if not required_cols.issubset(set(data_info.columns)):
        raise ValueError('plaque_sample_map entries must include sample_path/sample_id/genotype/region_location')

    data_info['key'] = data_info['sample_path'].apply(make_key)
    key_to_sid = dict(zip(data_info['key'], data_info['sample_id']))
    pdata.obs['sample_id'] = pdata.obs['key'].map(key_to_sid)
    pdata = pdata[pdata.obs['sample_id'].notna()].copy()

    if 'spatial' in pdata.obsm and pdata.obsm['spatial'].shape[1] >= 2:
        sp = np.asarray(pdata.obsm['spatial'])
        pdata.obsm['spatial_plot'] = sp[:, :2][:, [1, 0]]
    else:
        raise ValueError('Plaque spatial coordinates unavailable')

    # overlay
    ensure_parent(out_plaque_overlay)
    sids = data_info['sample_id'].astype(str).tolist()
    n = max(1, len(sids))
    ncols = 2
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5.5 * nrows))
    axes = np.asarray(axes).reshape(-1)
    for ax, sid in zip(axes, sids):
        ad = adata[adata.obs.get('sample_id', pd.Series(index=adata.obs_names, dtype=str)).astype(str) == sid]
        pd_sub = pdata[pdata.obs['sample_id'].astype(str) == sid]
        if ad.n_obs > 0 and 'center_x' in ad.obs.columns and 'center_y' in ad.obs.columns:
            ax.scatter(ad.obs['center_x'], ad.obs['center_y'], s=1, alpha=0.05, label='MERFISH')
        if pd_sub.n_obs > 0:
            vol = pd.to_numeric(pd_sub.obs.get('volume', pd.Series(np.ones(pd_sub.n_obs), index=pd_sub.obs_names)), errors='coerce').fillna(0).values
            if len(vol) > 0 and np.max(vol) > np.min(vol):
                size = ((vol - np.min(vol)) / (np.max(vol) - np.min(vol))) * 50.0
            else:
                size = np.ones(len(vol)) * 20.0
            ax.scatter(pd_sub.obsm['spatial_plot'][:, 0], pd_sub.obsm['spatial_plot'][:, 1], s=size, alpha=0.7, label='Plaque')
        info = data_info[data_info['sample_id'] == sid]
        if len(info) > 0:
            genotype = info['genotype'].iloc[0]
            region = info['region_location'].iloc[0]
            disease = 'WT' if genotype == 'WT' else 'AD'
            ax.set_title(f'{disease} | {region}', fontsize=14)
        else:
            ax.set_title(sid)
        ax.invert_yaxis()
        ax.set_aspect('equal')
        ax.legend(fontsize=8)
    for ax in axes[len(sids):]:
        ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_plaque_overlay, dpi=200)
    plt.close(fig)

    # burden bar
    ensure_parent(out_plaque_bar)
    rows = []
    for sid in sids:
        ad = adata[adata.obs.get('sample_id', pd.Series(index=adata.obs_names, dtype=str)).astype(str) == sid]
        pd_sub = pdata[pdata.obs['sample_id'].astype(str) == sid]
        volume_sum = float(pd.to_numeric(pd_sub.obs.get('volume', pd.Series(dtype=float)), errors='coerce').fillna(0).sum())
        cell_n = max(1, int(ad.n_obs))
        info = data_info[data_info['sample_id'] == sid]
        genotype = info['genotype'].iloc[0] if len(info) > 0 else 'NA'
        region = info['region_location'].iloc[0] if len(info) > 0 else 'NA'
        rows.append({'sample_id': sid, 'genotype': genotype, 'region_location': region, 'plaque_volume_per_cell': volume_sum / cell_n})
    burden = pd.DataFrame(rows)
    burden['status'] = burden['genotype'].replace({'5xFAD': 'AD'})
    burden['group'] = burden['status'].astype(str) + ' ' + burden['region_location'].astype(str)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.bar(burden['group'], burden['plaque_volume_per_cell'])
    ax.set_title('Plaque for each sample')
    ax.set_ylabel('Calculation of plaque')
    ax.tick_params(axis='x', rotation=45)
    fig.tight_layout()
    fig.savefig(out_plaque_bar, dpi=200)
    plt.close(fig)
except Exception as e:
    save_placeholder(out_plaque_overlay, f'Plaque overlay unavailable: {str(e)[:100]}')
    save_placeholder(out_plaque_bar, f'Plaque bar unavailable: {str(e)[:100]}')


# 6) Ligand-Receptor power vs Trem2/Apoe scatter (from DESeq2 VSD tables)
try:
    oligo = pd.read_csv(oligo_vsd_path, index_col=0)
    micro = pd.read_csv(micro_vsd_path, index_col=0)
    if not {'Ntm'}.issubset(set(oligo.index)) or not {'Stab1', 'Trem2', 'Apoe'}.issubset(set(micro.index)):
        raise ValueError('Required genes (Ntm, Stab1, Trem2, Apoe) missing in VSD tables')
    x = np.sqrt(pd.to_numeric(oligo.loc['Ntm'], errors='coerce') * pd.to_numeric(micro.loc['Stab1'], errors='coerce'))
    y = (pd.to_numeric(micro.loc['Trem2'], errors='coerce') + pd.to_numeric(micro.loc['Apoe'], errors='coerce')) / 2.0
    df = pd.DataFrame({'Ligand_Receptor_power': x, 'Trem2_Apoe_mean_expression': y}).dropna()
    if len(df) < 2:
        raise ValueError('Too few paired points for LR scatter')
    coef = np.polyfit(df['Ligand_Receptor_power'].values, df['Trem2_Apoe_mean_expression'].values, 1)
    xs = np.linspace(float(df['Ligand_Receptor_power'].min()), float(df['Ligand_Receptor_power'].max()), 100)
    ys = coef[0] * xs + coef[1]

    ensure_parent(out_lr_scatter)
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    ax.scatter(df['Ligand_Receptor_power'], df['Trem2_Apoe_mean_expression'], s=70)
    ax.plot(xs, ys, color='red', linestyle='--', linewidth=2)
    ax.set_xlabel('Ligand-Receptor power (sqrt(Ntm * Stab1))')
    ax.set_ylabel('(Trem2 + Apoe) / 2 on Microglia')
    ax.set_title('Ligand-Receptor vs Trem2/Apoe')
    fig.tight_layout()
    fig.savefig(out_lr_scatter, dpi=200)
    plt.close(fig)
except Exception as e:
    save_placeholder(out_lr_scatter, f'LR scatter unavailable: {str(e)[:100]}')


# 7) UMAP Stab1/Ntm and manual annotation
def find_gene_case_insensitive(var_names, target):
    tv = target.upper()
    for g in var_names:
        if str(g).upper() == tv:
            return g
    return None

try:
    if 'X_umap' not in adata.obsm:
        raise ValueError('X_umap not found in adata.obsm')
    g_stab1 = find_gene_case_insensitive(adata.var_names, 'Stab1')
    g_ntm = find_gene_case_insensitive(adata.var_names, 'Ntm')
    if g_stab1 is None or g_ntm is None:
        raise ValueError('Stab1/Ntm not found in var_names')
    sc.pl.umap(adata, color=[g_stab1, g_ntm], use_raw=False, show=False)
    ensure_parent(out_umap_stab1_ntm)
    plt.savefig(out_umap_stab1_ntm, dpi=200, bbox_inches='tight')
    plt.close('all')
except Exception as e:
    save_placeholder(out_umap_stab1_ntm, f'UMAP Stab1/Ntm unavailable: {str(e)[:100]}')

try:
    ann_col = pick_col(obs, ['cluster_pred', 'final_label_transfer', 'mannual_annotaion', 'Mannual_annotaion'])
    if 'X_umap' not in adata.obsm or ann_col is None:
        raise ValueError('X_umap or annotation column missing')
    sc.pl.umap(adata, color=[ann_col], use_raw=False, show=False)
    ensure_parent(out_umap_manual)
    plt.savefig(out_umap_manual, dpi=200, bbox_inches='tight')
    plt.close('all')
except Exception as e:
    save_placeholder(out_umap_manual, f'UMAP annotation unavailable: {str(e)[:100]}')
