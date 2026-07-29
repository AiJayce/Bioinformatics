configfile: 'config/config.yaml'

import os


def resolve_path(value):
    """Expand env/user markers and normalize to absolute path for source inputs."""
    if not isinstance(value, str):
        return value
    expanded = os.path.expanduser(os.path.expandvars(value))
    return os.path.abspath(expanded)


def resolve_dict(dct):
    return {k: resolve_path(v) for k, v in dct.items()}

ASSIGNMENT1 = config['assignment1']
SOURCE = resolve_dict(config['source'])
COPY = config['copy_targets']
SCVI_CFG = config.get('scvi', {})
RUNTIME_CFG = config.get('runtime', {})


def required_source_inputs():
    keys = [
        'raw_merscope_h5ad',
        'scvi_reference_h5ad',
        'latent_top100',
        'deseq2_dir',
        'nichenet_dir',
        'fgsea_dir',
        'vae_model_dir',
        'plaques_h5ad',
    ]
    return [SOURCE[k] for k in keys if SOURCE.get(k)]

rule all:
    input:
        'results/check_inputs.done',
        'results/scvi/annotated_scvi.h5ad',
        'results/scrna/preprocessed_scRNA.h5ad',
        'results/scrna/qc_violin(scRNA).png',
        'results/scrna/umap(scRNA).png',
        'results/scrna/leiden(scRNA).png',
        'results/assignment/after_bulk_ccf_annotated.h5ad',
        'results/assignment/after_bulk_ccf_obs.csv',
        'results/assignment/basic_qc_violin(scRNA).png',
        'results/assignment/ccf_top20(scRNA).png',
        'results/assignment/bulk_ccf_wt_vs_ad(scRNA).png',
        'results/assignment/scanpy_umap_summary(scRNA).png',
        'results/assignment/scanpy_pca_summary(scRNA).png',
        'results/assignment/cluster_plot(scRNA).png',
        'results/assignment/leiden_plot(scRNA).png',
        'results/assignment/notebook_pca_elbow(scRNA).png',
        'results/assignment/notebook_doublet_score_hist(scRNA).png',
        'results/assignment/notebook_scvi_loss_curve(scRNA).png',
        'results/assignment/notebook_rank_genes_dotplot(scRNA).png',
        'results/assignment/notebook_ccf_registration_qc(scRNA).png',
        'results/assignment/notebook_bbknn_harmony_umap_compare(scRNA).png',
        'results/assignment/notebook_class_score_violin(scRNA).png',
        COPY['latent_top100_csv'],
        *COPY['deseq2_outputs'],
        'results/assignment2/latent_dim_gene_count(scRNA).png',
        'results/assignment2/latent_effect_distribution(scRNA).png',
        'results/assignment2/latent_top_genes_heatmap(scRNA).png',
        'results/assignment2/z15_top_genes_bar(scRNA).png',
        'results/assignment2/microglia_vsd_pca.png',
        'results/assignment2/oligodendrocyte_vsd_pca.png',
        'results/assignment2/microglia_neighbor_assoc_heatmap(scRNA).png',
        'results/assignment2/microglia_trem2_vs_oligo_neighbors(scRNA).png',
        'results/assignment2/microglia_apoe_vs_oligo_neighbors(scRNA).png',
        'results/assignment2/microglia_pseudobulk_heatmap.png',
        'results/assignment2/oligodendrocyte_pseudobulk_heatmap.png',
        'results/assignment2/notebook_trem2_apoe_boxplot(scRNA).png',
        'results/assignment2/notebook_volcano_wt_posterior_vs_anterior(scRNA).png',
        'results/assignment2/notebook_volcano_ad_posterior_vs_anterior(scRNA).png',
        'results/assignment2/notebook_venn_deg_up(scRNA).png',
        'results/assignment2/notebook_spatial_apoe_trem2_overlay(scRNA).png',
        'results/assignment2/notebook_spatial_apoe_trem2_zoom(scRNA).png',
        'results/assignment2/notebook_vae_latent_heatmap(scRNA).png',
        'results/assignment2/notebook_vae_z15_top10(scRNA).png',
        'results/assignment2/notebook_vae_top10_distribution(scRNA).png',
        'results/assignment2/notebook_plaque_merfish_overlay(scRNA).png',
        'results/assignment2/notebook_plaque_per_sample_bar(scRNA).png',
        'results/assignment2/notebook_lr_power_vs_trem2_apoe(scRNA).png',
        'results/assignment2/notebook_umap_stab1_ntm(scRNA).png',
        'results/assignment2/notebook_umap_manual_annotation(scRNA).png',
        *COPY['nichenet_outputs'],
        COPY['fgsea_outputs']['results_csv'],
        COPY['fgsea_outputs']['plot_png'],
        'results/external/deseq2_sync.done',
        'results/external/nichenet_sync.done',
        'results/external/fgsea_sync.done',
        'results/combined/analysis_complete.txt'

rule scrna_preprocess:
    input:
        'results/check_inputs.done',
        source='results/scvi/annotated_scvi.h5ad',
        config='config/config.yaml'
    output:
        h5ad='results/scrna/preprocessed_scRNA.h5ad',
        qc='results/scrna/qc_violin(scRNA).png',
        umap='results/scrna/umap(scRNA).png',
        leiden='results/scrna/leiden(scRNA).png'
    script:
        'scripts/00b_scrna_preprocess.py'

rule scvi_annotate_raw:
    input:
        'results/check_inputs.done',
        query_raw=SOURCE['raw_merscope_h5ad'],
        reference=SOURCE['scvi_reference_h5ad']
    output:
        annotated='results/scvi/annotated_scvi.h5ad',
        history='results/scvi/scvi_history.json'
    params:
        scvi_python=resolve_path(RUNTIME_CFG.get('scvi_python', 'python')),
        label_key=SCVI_CFG.get('label_key', 'final_label_transfer'),
        batch_key=SCVI_CFG.get('batch_key', 'sample_id'),
        unlabeled_category=SCVI_CFG.get('unlabeled_category', 'Unknown'),
        n_latent=SCVI_CFG.get('n_latent', 30),
        max_epochs_scvi=SCVI_CFG.get('max_epochs_scvi', 200),
        max_epochs_scanvi=SCVI_CFG.get('max_epochs_scanvi', 40),
        seed=SCVI_CFG.get('seed', 0)
    shell:
        (
            "{params.scvi_python} scripts/00a_scvi_annotate.py "
            "--query-h5ad '{input.query_raw}' "
            "--reference-h5ad '{input.reference}' "
            "--output-h5ad '{output.annotated}' "
            "--history-json '{output.history}' "
            "--label-key '{params.label_key}' "
            "--batch-key '{params.batch_key}' "
            "--unlabeled-category '{params.unlabeled_category}' "
            "--n-latent {params.n_latent} "
            "--max-epochs-scvi {params.max_epochs_scvi} "
            "--max-epochs-scanvi {params.max_epochs_scanvi} "
            "--seed {params.seed}"
        )

rule check_inputs:
    input:
        required=required_source_inputs()
    output:
        'results/check_inputs.done'
    script:
        'scripts/00_check_inputs.py'

rule assignment1_ccf:
    input:
        'results/check_inputs.done',
        h5ad='results/scrna/preprocessed_scRNA.h5ad',
        config='config/config.yaml'
    output:
        h5ad='results/assignment/after_bulk_ccf_annotated.h5ad',
        csv='results/assignment/after_bulk_ccf_obs.csv',
        basic_qc='results/assignment/basic_qc_violin(scRNA).png',
        ccf_top20='results/assignment/ccf_top20(scRNA).png',
        bulk_png='results/assignment/bulk_ccf_wt_vs_ad(scRNA).png',
        scanpy_umap='results/assignment/scanpy_umap_summary(scRNA).png',
        scanpy_pca='results/assignment/scanpy_pca_summary(scRNA).png',
        cluster_plot='results/assignment/cluster_plot(scRNA).png',
        leiden_plot='results/assignment/leiden_plot(scRNA).png',
        summary='results/assignment/pipeline_summary.txt'
    script:
        'scripts/01_assignment1_real.py'

rule assignment1_notebook_figures:
    input:
        adata='results/assignment/after_bulk_ccf_annotated.h5ad',
        config='config/config.yaml'
    output:
        pca_elbow='results/assignment/notebook_pca_elbow(scRNA).png',
        doublet_hist='results/assignment/notebook_doublet_score_hist(scRNA).png',
        scvi_loss='results/assignment/notebook_scvi_loss_curve(scRNA).png',
        rank_dotplot='results/assignment/notebook_rank_genes_dotplot(scRNA).png',
        ccf_qc='results/assignment/notebook_ccf_registration_qc(scRNA).png',
        bbknn_harmony='results/assignment/notebook_bbknn_harmony_umap_compare(scRNA).png',
        class_score_violin='results/assignment/notebook_class_score_violin(scRNA).png'
    script:
        'scripts/08_assignment1_notebook_figures.py'

rule assignment2_latent:
    input:
        'results/check_inputs.done',
        assignment1_h5ad='results/assignment/after_bulk_ccf_annotated.h5ad',
        source=SOURCE['latent_top100']
    output:
        COPY['latent_top100_csv']
    script:
        'scripts/02_copy_file.py'

rule deseq2_collect:
    input:
        'results/check_inputs.done',
        assignment1_h5ad='results/assignment/after_bulk_ccf_annotated.h5ad',
        source_dir=SOURCE['deseq2_dir']
    output:
        COPY['deseq2_outputs']
    script:
        'scripts/03_copy_many_files.py'

rule nichenet_collect:
    input:
        'results/check_inputs.done',
        assignment1_h5ad='results/assignment/after_bulk_ccf_annotated.h5ad',
        deseq2_vsd_oligo='results/deseq2/Oligodendrocyte_bulk_rna_vsd_preprocess.csv',
        deseq2_vsd_micro='results/deseq2/Microglia_bulk_rna_vsd_preprocess.csv',
        source_dir=SOURCE['nichenet_dir']
    output:
        COPY['nichenet_outputs']
    script:
        'scripts/03_copy_many_files.py'

rule assignment2_plots:
    input:
        latent=COPY['latent_top100_csv'],
        adata='results/assignment/after_bulk_ccf_annotated.h5ad',
        micro_vsd='results/deseq2/Microglia_bulk_rna_vsd_preprocess.csv',
        oligo_vsd='results/deseq2/Oligodendrocyte_bulk_rna_vsd_preprocess.csv',
        micro_pseudobulk='results/deseq2/pseudobulk_Microglia.csv',
        oligo_pseudobulk='results/deseq2/pseudobulk_Oligodendrocyte.csv'
    output:
        latent_dim_count='results/assignment2/latent_dim_gene_count(scRNA).png',
        effect_dist='results/assignment2/latent_effect_distribution(scRNA).png',
        latent_heatmap='results/assignment2/latent_top_genes_heatmap(scRNA).png',
        z15_bar='results/assignment2/z15_top_genes_bar(scRNA).png',
        micro_pca='results/assignment2/microglia_vsd_pca.png',
        oligo_pca='results/assignment2/oligodendrocyte_vsd_pca.png',
        neighbor_heatmap='results/assignment2/microglia_neighbor_assoc_heatmap(scRNA).png',
        trem2_scatter='results/assignment2/microglia_trem2_vs_oligo_neighbors(scRNA).png',
        apoe_scatter='results/assignment2/microglia_apoe_vs_oligo_neighbors(scRNA).png',
        micro_pseudobulk_heatmap='results/assignment2/microglia_pseudobulk_heatmap.png',
        oligo_pseudobulk_heatmap='results/assignment2/oligodendrocyte_pseudobulk_heatmap.png'
    script:
        'scripts/05_assignment2_plots.py'

rule assignment2_notebook_figures:
    input:
        adata='results/assignment/after_bulk_ccf_annotated.h5ad',
        latent=COPY['latent_top100_csv'],
        vae_model_dir=SOURCE['vae_model_dir'],
        plaques_h5ad=SOURCE['plaques_h5ad'],
        oligo_vsd='results/deseq2/Oligodendrocyte_bulk_rna_vsd_preprocess.csv',
        micro_vsd='results/deseq2/Microglia_bulk_rna_vsd_preprocess.csv',
        config='config/config.yaml'
    output:
        boxplot='results/assignment2/notebook_trem2_apoe_boxplot(scRNA).png',
        volcano_wt='results/assignment2/notebook_volcano_wt_posterior_vs_anterior(scRNA).png',
        volcano_ad='results/assignment2/notebook_volcano_ad_posterior_vs_anterior(scRNA).png',
        venn='results/assignment2/notebook_venn_deg_up(scRNA).png',
        spatial='results/assignment2/notebook_spatial_apoe_trem2_overlay(scRNA).png',
        spatial_zoom='results/assignment2/notebook_spatial_apoe_trem2_zoom(scRNA).png',
        vae_heatmap='results/assignment2/notebook_vae_latent_heatmap(scRNA).png',
        vae_top10='results/assignment2/notebook_vae_z15_top10(scRNA).png',
        vae_distribution='results/assignment2/notebook_vae_top10_distribution(scRNA).png',
        plaque_overlay='results/assignment2/notebook_plaque_merfish_overlay(scRNA).png',
        plaque_bar='results/assignment2/notebook_plaque_per_sample_bar(scRNA).png',
        lr_scatter='results/assignment2/notebook_lr_power_vs_trem2_apoe(scRNA).png',
        umap_stab1_ntm='results/assignment2/notebook_umap_stab1_ntm(scRNA).png',
        umap_manual='results/assignment2/notebook_umap_manual_annotation(scRNA).png'
    script:
        'scripts/06_assignment2_notebook_figures.py'

rule sync_deseq2_all:
    input:
        source_dir=SOURCE['deseq2_dir']
    output:
        done='results/external/deseq2_sync.done'
    params:
        name='DESeq2'
    script:
        'scripts/07_sync_external_dir.py'

rule sync_nichenet_all:
    input:
        source_dir=SOURCE['nichenet_dir']
    output:
        done='results/external/nichenet_sync.done'
    params:
        name='NicheNet'
    script:
        'scripts/07_sync_external_dir.py'

rule sync_fgsea_all:
    input:
        source_dir=SOURCE['fgsea_dir']
    output:
        done='results/external/fgsea_sync.done'
    params:
        name='fGSEA'
    script:
        'scripts/07_sync_external_dir.py'

rule fgsea_real:
    input:
        'results/check_inputs.done',
        assignment1_h5ad='results/assignment/after_bulk_ccf_annotated.h5ad',
        latent=COPY['latent_top100_csv']
    output:
        results_csv=COPY['fgsea_outputs']['results_csv'],
        plot_png=COPY['fgsea_outputs']['plot_png']
    shell:
        "Rscript 'scripts/04_fgsea_real.R' '{input.latent}' '{output.results_csv}' '{output.plot_png}'"

rule combined_complete:
    input:
        'results/check_inputs.done',
        'results/scrna/preprocessed_scRNA.h5ad',
        'results/scrna/qc_violin(scRNA).png',
        'results/scrna/umap(scRNA).png',
        'results/scrna/leiden(scRNA).png',
        'results/assignment/after_bulk_ccf_annotated.h5ad',
        'results/assignment/after_bulk_ccf_obs.csv',
        'results/assignment/basic_qc_violin(scRNA).png',
        'results/assignment/ccf_top20(scRNA).png',
        'results/assignment/bulk_ccf_wt_vs_ad(scRNA).png',
        'results/assignment/scanpy_umap_summary(scRNA).png',
        'results/assignment/scanpy_pca_summary(scRNA).png',
        'results/assignment/cluster_plot(scRNA).png',
        'results/assignment/leiden_plot(scRNA).png',
        'results/assignment/notebook_pca_elbow(scRNA).png',
        'results/assignment/notebook_doublet_score_hist(scRNA).png',
        'results/assignment/notebook_scvi_loss_curve(scRNA).png',
        'results/assignment/notebook_rank_genes_dotplot(scRNA).png',
        'results/assignment/notebook_ccf_registration_qc(scRNA).png',
        'results/assignment/notebook_bbknn_harmony_umap_compare(scRNA).png',
        'results/assignment/notebook_class_score_violin(scRNA).png',
        COPY['latent_top100_csv'],
        *COPY['deseq2_outputs'],
        'results/assignment2/latent_dim_gene_count(scRNA).png',
        'results/assignment2/latent_effect_distribution(scRNA).png',
        'results/assignment2/latent_top_genes_heatmap(scRNA).png',
        'results/assignment2/z15_top_genes_bar(scRNA).png',
        'results/assignment2/microglia_vsd_pca.png',
        'results/assignment2/oligodendrocyte_vsd_pca.png',
        'results/assignment2/microglia_neighbor_assoc_heatmap(scRNA).png',
        'results/assignment2/microglia_trem2_vs_oligo_neighbors(scRNA).png',
        'results/assignment2/microglia_apoe_vs_oligo_neighbors(scRNA).png',
        'results/assignment2/microglia_pseudobulk_heatmap.png',
        'results/assignment2/oligodendrocyte_pseudobulk_heatmap.png',
        'results/assignment2/notebook_trem2_apoe_boxplot(scRNA).png',
        'results/assignment2/notebook_volcano_wt_posterior_vs_anterior(scRNA).png',
        'results/assignment2/notebook_volcano_ad_posterior_vs_anterior(scRNA).png',
        'results/assignment2/notebook_venn_deg_up(scRNA).png',
        'results/assignment2/notebook_spatial_apoe_trem2_overlay(scRNA).png',
        'results/assignment2/notebook_spatial_apoe_trem2_zoom(scRNA).png',
        'results/assignment2/notebook_vae_latent_heatmap(scRNA).png',
        'results/assignment2/notebook_vae_z15_top10(scRNA).png',
        'results/assignment2/notebook_vae_top10_distribution(scRNA).png',
        'results/assignment2/notebook_plaque_merfish_overlay(scRNA).png',
        'results/assignment2/notebook_plaque_per_sample_bar(scRNA).png',
        'results/assignment2/notebook_lr_power_vs_trem2_apoe(scRNA).png',
        'results/assignment2/notebook_umap_stab1_ntm(scRNA).png',
        'results/assignment2/notebook_umap_manual_annotation(scRNA).png',
        *COPY['nichenet_outputs'],
        COPY['fgsea_outputs']['results_csv'],
        COPY['fgsea_outputs']['plot_png'],
        'results/external/deseq2_sync.done',
        'results/external/nichenet_sync.done',
        'results/external/fgsea_sync.done'
    output:
        'results/combined/analysis_complete.txt'
    shell:
        "mkdir -p results/combined && printf 'Combined workflow complete\n' > {output}"
