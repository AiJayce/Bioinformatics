args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop('Usage: Rscript 04_fgsea_real.R <latent_csv> <results_csv> <plot_png>')
}

latent_csv <- args[[1]]
results_csv <- args[[2]]
plot_png <- args[[3]]

suppressPackageStartupMessages({
  library(ggplot2)
  library(msigdbr)
  library(fgsea)
  library(stringr)
})

df <- read.csv(latent_csv, stringsAsFactors = FALSE)

if (!'gene' %in% names(df) && ncol(df) >= 1) {
  names(df)[1] <- 'gene'
}
if (!'latent_dim' %in% names(df) && ncol(df) >= 2) {
  names(df)[2] <- 'latent_dim'
}
if (!'effect' %in% names(df) && ncol(df) >= 3) {
  names(df)[3] <- 'effect'
}

df$gene <- as.character(df$gene)
df$latent_dim <- as.character(df$latent_dim)
df$effect <- as.numeric(df$effect)

mouse_go_bp <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'GO:BP')
mouse_go_cc <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'GO:CC')
mouse_go_mf <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'GO:MF')
mouse_reactome <- msigdbr(species = 'Mus musculus', category = 'C2', subcategory = 'CP:REACTOME')
mouse_hallmark <- msigdbr(species = 'Mus musculus', category = 'H')

mouse_pathway <- rbind(mouse_go_bp, mouse_go_cc, mouse_go_mf, mouse_reactome, mouse_hallmark)
gene_sets <- split(mouse_pathway$gene_symbol, paste(mouse_pathway$gs_cat, mouse_pathway$gs_name, sep = '_'))

normalize_fgsea_df <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  expected <- c('pathway', 'pval', 'padj', 'log2err', 'ES', 'NES', 'size', 'leadingEdge')
  if (ncol(x) == 0) {
    out <- as.data.frame(setNames(replicate(length(expected), vector(), simplify = FALSE), expected), stringsAsFactors = FALSE)
    return(out)
  }
  n <- min(length(expected), ncol(x))
  names(x)[seq_len(n)] <- expected[seq_len(n)]
  x
}

empty_fgsea_df <- data.frame(
  pathway = character(),
  pval = numeric(),
  padj = numeric(),
  log2err = numeric(),
  ES = numeric(),
  NES = numeric(),
  size = numeric(),
  leadingEdge = I(list()),
  latent_dim = character(),
  stringsAsFactors = FALSE
)

fgsea_results <- list()
for (ld in unique(df$latent_dim)) {
  message('Running ', ld)
  stats_df <- df[df$latent_dim == ld, c('gene', 'effect')]
  stats_df <- stats_df[!duplicated(stats_df$gene), , drop = FALSE]
  stats <- setNames(stats_df$effect, stats_df$gene)
  stats <- sort(stats, decreasing = TRUE)
  res <- fgsea(pathways = gene_sets, stats = stats, minSize = 10, maxSize = 500)
  fgsea_results[[ld]] <- res
}

fgsea_results_df <- do.call(rbind, lapply(names(fgsea_results), function(ld) {
  x <- normalize_fgsea_df(fgsea_results[[ld]])
  x$latent_dim <- ld
  x
}))
fgsea_results_df <- as.data.frame(fgsea_results_df, stringsAsFactors = FALSE)
if (ncol(fgsea_results_df) == 0) {
  fgsea_results_df <- empty_fgsea_df
}
if ('leadingEdge' %in% names(fgsea_results_df) && is.list(fgsea_results_df$leadingEdge)) {
  fgsea_results_df$leadingEdge <- vapply(fgsea_results_df$leadingEdge, function(x) {
    if (length(x) == 0) {
      ''
    } else {
      paste(x, collapse = ';')
    }
  }, character(1))
}
write.csv(fgsea_results_df, results_csv, row.names = FALSE)

target_pathway <- c(
  'C5_GOBP_INFLAMMATORY_RESPONSE',
  'C5_GOBP_REGULATION_OF_CELL_ACTIVATION',
  'C5_GOBP_REGULATION_OF_LIPID_METABOLIC_PROCESS',
  'C2_REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM',
  'C5_GOBP_CELL_ACTIVATION'
)

plot_df <- if ('z15' %in% names(fgsea_results)) {
  x <- normalize_fgsea_df(fgsea_results[['z15']])
  x$latent_dim <- 'z15'
  if (ncol(x) == 0) {
    empty_fgsea_df
  } else {
    names(x)[1] <- 'pathway'
    x[x[[1]] %in% target_pathway, , drop = FALSE]
  }
} else {
  if (ncol(fgsea_results_df) == 0) {
    empty_fgsea_df
  } else {
    names(fgsea_results_df)[1] <- 'pathway'
    fgsea_results_df[fgsea_results_df[[1]] %in% target_pathway, , drop = FALSE]
  }
}

if (nrow(plot_df) == 0) {
  if (ncol(fgsea_results_df) == 0) {
    plot_df <- empty_fgsea_df
  } else {
    names(fgsea_results_df)[1] <- 'pathway'
    plot_df <- fgsea_results_df[fgsea_results_df[[1]] %in% target_pathway, , drop = FALSE]
  }
}

plot_df$overlap_percent <- mapply(function(p, s) { s / length(gene_sets[[p]]) * 100 }, plot_df$pathway, plot_df$size)
plot_df$pathway_label <- str_to_title(str_replace_all(str_replace(str_replace(plot_df$pathway, '^C5_GOBP_', ''), '^C2_REACTOME_', ''), '_', ' '))
plot_df$significance <- ifelse(plot_df$pval < 0.05, '*', 'n.s.')
plot_df$neg_log10_p <- -log10(plot_df$pval)
plot_df <- plot_df[order(-plot_df$NES), , drop = FALSE]

plot_df$pathway_label <- factor(plot_df$pathway_label, levels = rev(plot_df$pathway_label))

png(plot_png, width = 1800, height = 1000, res = 180, type = 'cairo')
p <- ggplot(plot_df, aes(x = NES, y = pathway_label, size = overlap_percent)) +
  geom_point(aes(color = neg_log10_p)) +
  scale_color_gradient(low = 'lightgrey', high = 'blue', name = '-log10(p-value)') +
  geom_text(aes(label = significance), nudge_y = 0.35, color = 'red', fontface = 'bold', size = 6, show.legend = FALSE) +
  scale_size(range = c(8, 12), name = 'Gene overlap (%)') +
  theme_classic() +
  labs(x = 'Enrichment Score (NES)', y = NULL)
print(p)
dev.off()
