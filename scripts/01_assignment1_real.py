import os
import yaml
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import nrrd


def resolve_path(path_str):
    if not isinstance(path_str, str):
        return path_str
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path_str)))


def load_micron2px_matrix(csv_path):
    m = np.loadtxt(csv_path)
    if m.shape != (3, 3):
        raise ValueError(f'Expected 3x3 matrix, got {m.shape} from {csv_path}')
    return m


def apply_affine_2d(x, y, m3x3):
    pts = np.column_stack([x, y, np.ones(len(x))])
    out = pts @ m3x3.T
    return out[:, 0], out[:, 1]


def apply_2x3(x, y, m2x3):
    pts = np.column_stack([x, y, np.ones(len(x))])
    out = pts @ np.asarray(m2x3, dtype=float).T
    return out[:, 0], out[:, 1]


cfg = getattr(snakemake, 'config', None)
if not isinstance(cfg, dict) or len(cfg) == 0:
    with open(snakemake.input.config, 'r', encoding='utf-8') as f:
        cfg = yaml.safe_load(f)

cfg1 = cfg['assignment1']
data_dir = resolve_path(cfg['source']['atlas_data_dir'])
input_h5ad = snakemake.input.h5ad if hasattr(snakemake.input, 'h5ad') else resolve_path(cfg1['existing_h5ad_path'])
adata = sc.read_h5ad(input_h5ad)

annotation_vol, _ = nrrd.read(resolve_path(cfg1['atlas_annot_path']))
structure_tree = pd.read_csv(resolve_path(cfg1['structure_tree_path']))

id_col = 'id' if 'id' in structure_tree.columns else structure_tree.columns[0]
acr_col = 'acronym' if 'acronym' in structure_tree.columns else None
name_col = 'name' if 'name' in structure_tree.columns else None

id_to_acr = {}
id_to_name = {}
for _, r in structure_tree.iterrows():
    rid = int(r[id_col])
    id_to_acr[rid] = str(r[acr_col]) if acr_col is not None else str(rid)
    id_to_name[rid] = str(r[name_col]) if name_col is not None else str(rid)

adata.obs['x_mosaic_px'] = np.nan
adata.obs['y_mosaic_px'] = np.nan
adata.obs['x_ccf_vox'] = np.nan
adata.obs['y_ccf_vox'] = np.nan
adata.obs['z_ccf_vox'] = np.nan
adata.obs['ccf_reg_mode'] = 'none'

for sid in adata.obs['sample_id'].astype(str).unique():
    smask = adata.obs['sample_id'].astype(str) == sid
    t_csv = os.path.join(data_dir, cfg1['sample_to_micron2px_csv'][sid])
    m_micron2px = load_micron2px_matrix(t_csv)
    x_um = adata.obs.loc[smask, 'center_x'].values
    y_um = adata.obs.loc[smask, 'center_y'].values
    x_px, y_px = apply_affine_2d(x_um, y_um, m_micron2px)
    adata.obs.loc[smask, 'x_mosaic_px'] = x_px
    adata.obs.loc[smask, 'y_mosaic_px'] = y_px
    reg = cfg1['sample_registration'][sid]
    x_vox, y_vox = apply_2x3(x_px, y_px, reg['matrix_2x3'])
    adata.obs.loc[smask, 'x_ccf_vox'] = x_vox
    adata.obs.loc[smask, 'y_ccf_vox'] = y_vox
    adata.obs.loc[smask, 'z_ccf_vox'] = float(reg['z_voxel'])
    adata.obs.loc[smask, 'ccf_reg_mode'] = 'inferred_coronal_2d'

z_max, y_max, x_max = annotation_vol.shape
adata.obs['ccf_id'] = -1
adata.obs['CCF_acronym'] = 'Not_mapped'
adata.obs['CCF_name'] = 'Not_mapped'

valid = adata.obs[['x_ccf_vox', 'y_ccf_vox', 'z_ccf_vox']].notna().all(axis=1)
xv = np.rint(adata.obs.loc[valid, 'x_ccf_vox'].values).astype(int)
yv = np.rint(adata.obs.loc[valid, 'y_ccf_vox'].values).astype(int)
zv = np.rint(adata.obs.loc[valid, 'z_ccf_vox'].values).astype(int)
in_bounds = (xv >= 0) & (xv < x_max) & (yv >= 0) & (yv < y_max) & (zv >= 0) & (zv < z_max)
idx_in = adata.obs.index[valid][in_bounds]
ann = annotation_vol[zv[in_bounds], yv[in_bounds], xv[in_bounds]].astype(int)
adata.obs.loc[idx_in, 'ccf_id'] = ann
adata.obs.loc[idx_in, 'CCF_acronym'] = [id_to_acr.get(int(i), 'Unknown') for i in ann]
adata.obs.loc[idx_in, 'CCF_name'] = [id_to_name.get(int(i), 'Unknown') for i in ann]

qc_cols = ['total_counts', 'n_genes_by_counts', 'pct_counts_mt']
if all(c in adata.obs.columns for c in ['total_counts', 'n_genes_by_counts']):
    if 'mt' not in adata.var.columns:
        adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    if 'ribo' not in adata.var.columns:
        adata.var['ribo'] = adata.var_names.str.upper().str.startswith(('RPS', 'RPL'))
    if 'hb' not in adata.var.columns:
        adata.var['hb'] = adata.var_names.str.upper().str.startswith(('HBA', 'HBB'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], percent_top=None, inplace=True)

fig, ax = plt.subplots(figsize=(7, 4))
counts = adata.obs['CCF_acronym'].value_counts().head(20)
counts.plot(kind='bar', ax=ax)
ax.set_title('Top 20 CCF labels')
ax.set_xlabel('CCF_acronym')
ax.set_ylabel('count')
plt.xticks(rotation=70)
plt.tight_layout()
fig.savefig(snakemake.output.ccf_top20, dpi=200)
plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 4))
ax.scatter(adata.obs['center_x'], adata.obs['center_y'], s=1, alpha=0.2)
ax.set_title('Basic spatial QC')
ax.set_xlabel('center_x')
ax.set_ylabel('center_y')
ax.invert_yaxis()
plt.tight_layout()
fig.savefig(snakemake.output.basic_qc, dpi=200)
plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 4))
wt = adata[adata.obs['sample_id'].astype(str).str.contains('WT')]
ad = adata[~adata.obs['sample_id'].astype(str).str.contains('WT')]
wt_prop = (wt.obs['CCF_acronym'].value_counts(normalize=True).head(10))
ad_prop = (ad.obs['CCF_acronym'].value_counts(normalize=True).head(10))
plot_df = pd.DataFrame({'WT': wt_prop, 'AD': ad_prop}).fillna(0)
plot_df.plot(kind='barh', ax=ax)
ax.set_title('Bulk CCF WT vs AD')
plt.tight_layout()
fig.savefig(snakemake.output.bulk_png, dpi=200)
plt.close(fig)

def save_scanpy_embedding(kind, out_path, color_keys, fallback_title):
    has_embedding = (f'X_{kind}' in adata.obsm) if hasattr(adata, 'obsm') else False
    keys = [c for c in color_keys if c in adata.obs.columns]
    if has_embedding and keys:
        getattr(sc.pl, kind)(adata, color=keys, ncols=2, show=False)
        plt.savefig(out_path, dpi=200, bbox_inches='tight')
        plt.close('all')
    else:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, fallback_title, ha='center', va='center')
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(out_path, dpi=200)
        plt.close(fig)

save_scanpy_embedding('umap', snakemake.output.scanpy_umap, ['sample_id', 'CCF_acronym', 'CCF_major', 'final_label_transfer'], 'UMAP unavailable')
save_scanpy_embedding('pca', snakemake.output.scanpy_pca, ['sample_id', 'CCF_major'], 'PCA unavailable')

cluster_keys = [c for c in adata.obs.columns if c.startswith('cluster') or c == 'final_label_transfer']
if 'X_umap' in adata.obsm and cluster_keys:
    sc.pl.umap(adata, color=cluster_keys[:4], ncols=2, show=False)
    plt.savefig(snakemake.output.cluster_plot, dpi=200, bbox_inches='tight')
    plt.close('all')
else:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, 'Cluster plot unavailable', ha='center', va='center')
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(snakemake.output.cluster_plot, dpi=200)
    plt.close(fig)

leiden_keys = [c for c in adata.obs.columns if c.startswith('leiden_res_')]
if 'X_umap' in adata.obsm and leiden_keys:
    sc.pl.umap(adata, color=leiden_keys[:4], ncols=2, show=False)
    plt.savefig(snakemake.output.leiden_plot, dpi=200, bbox_inches='tight')
    plt.close('all')
else:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, 'Leiden plot unavailable', ha='center', va='center')
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(snakemake.output.leiden_plot, dpi=200)
    plt.close(fig)

adata.obs.to_csv(snakemake.output.csv)
adata.write_h5ad(snakemake.output.h5ad)

with open(snakemake.output.summary, 'w', encoding='utf-8') as f:
    f.write(f'n_obs: {adata.n_obs}\n')
    f.write(f'mapped: {(adata.obs["CCF_acronym"].astype(str) != "Not_mapped").sum()}\n')
