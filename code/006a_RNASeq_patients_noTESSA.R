library(EnhancedVolcano)
library(RColorBrewer)
library(ComplexHeatmap)
library(DESeq2)
library(msigdbr)
library(patchwork)
library(tidyverse)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")

# DESeq2 Differential Expression ----
if(!all(file.exists("data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_res_shrunk.rds", 
                    "data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_dds.rds"))){
  expmat_orig <- read_table("data/RNASeq_Lasorsa/BGI_F21FTSEUHT2056_HOMypeiR/Expression/gene_expression.xls") %>% 
    select(gene_id, gene_symbol, starts_with("read_count"), -contains(c("KO52", "ctrTESSA"))) 
  
  expmat <- expmat_orig %>% 
    select(-gene_id) %>% 
    column_to_rownames("gene_symbol") %>% 
    rename_with(~gsub("read_count_", "", .x, fixed = TRUE)) %>% 
    mutate(across(everything(), ~round(.x, 0))) %>% 
    as.matrix()
  
  coldata <- tibble(sample = colnames(expmat),
                    genotype = ifelse(str_starts(sample, "agc"), "kd", "wt"),
                    cell_line = str_remove(sample, "_[1,2,3]$")) %>% 
    column_to_rownames("sample") %>% 
    as.matrix()
  
  all(rownames(coldata) %in% colnames(expmat))
  all(rownames(coldata) == colnames(expmat))
  
  dds <- DESeqDataSetFromMatrix(countData = expmat,
                                colData = coldata,
                                design = ~genotype)
  
  dds$genotype <- relevel(dds$genotype, ref = "wt")
  
  # Dimensionality reduction and sample correlation
  rld <- rlog(dds, blind = FALSE)
  sampleDists <- dist(t(assay(rld)))    # Euclidean distance between samples
  sampleDistMatrix <- as.matrix(sampleDists)
  colnames(sampleDistMatrix) <- NULL
  rownames(sampleDistMatrix) <- paste0(rep(c("P1_", "P2A_", "P2B_", "C1_", "C2_"), each = 3), 1:3)
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  png("plots/006a_RNASeq_patients_noTESSA/006a_1_Euclidean_distance.png", h = 2200, w = 3500, res = 600)
  Heatmap(sampleDistMatrix,
          col = colors,
          column_title = "Overall similarity between samples",
          name = "Euclidean \ndistance\n",
          rect_gp = gpar(col = "grey60", lwd = 1)
  )
  dev.off()
  
  plotPCA(rld, intgroup = c("cell_line"), ntop = nrow(rld)) +
    geom_point(aes(fill = group), size = 5, shape = 21, color = "black") +
    labs(title = "Principal Component Analysis") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    guides(color = "none") +
    scale_fill_hue(labels = c("P1", "P2A", "P2B", "C1", "C2"))
  
  ggsave("plots/006a_RNASeq_patients_noTESSA/006a_2_PCA_plot.png",
         h = 1500,
         w = 2000,
         units = "px",
         dpi = 300)
  
  # Repeat graph but with means of samples
  rld_mean <- matrix(nrow = nrow(assay(rld)),
                     ncol = 5,
                     data = NA,
                     dimnames = list(rownames(assay(rld)), c("P1", "P2A", "P2B", "C1", "C2")))
  rld_mean[,"P1"] <- rowMeans(assay(rld)[,1:3])
  rld_mean[,"P2A"] <- rowMeans(assay(rld)[,4:6])
  rld_mean[,"P2B"] <- rowMeans(assay(rld)[,7:9])
  rld_mean[,"C1"] <- rowMeans(assay(rld)[,10:12])
  rld_mean[,"C2"] <- rowMeans(assay(rld)[,13:15])
  
  hm_mean <- dist(t(rld_mean)) |> 
    as.matrix()
  
  png("plots/006a_RNASeq_patients_noTESSA/006_1_Euclidean_distance_mean.png", h = 1200, w = 2000, res = 400)
  Heatmap(hm_mean,
          col = colors,
          column_title = "Overall similarity between samples",
          name = "Euclidean \ndistance\n",
          rect_gp = gpar(col = "grey60", lwd = 1)
  )
  dev.off()
  
  
  pca_rld_mean <- prcomp(t(rld_mean[which(rowSums(rld_mean) != 0),]), scale. = T)
  summary_pca <- summary(pca_rld_mean)
  pca_rld_mean$x |> 
    as_tibble(rownames = "Sample") |> 
    ggplot(aes(PC1, PC2, color = Sample)) +
    geom_point(size = 5) + 
    labs(title = "Principal Component Analysis",
         x = paste0("PC1 (", round(summary_pca$importance[2, "PC1"]*100, 2), "%)"),
         y = paste0("PC2 (", round(summary_pca$importance[2, "PC2"]*100, 2), "%)")) +
    theme_light()
  ggsave("plots/006a_RNASeq_patients_noTESSA/006a_2_PCA_plot_mean.png",
         h = 1000,
         w = 1500,
         units = "px",
         dpi = 300)
  
  # Differential expression
  des <- DESeq(dds)
  res <- results(des)
  resLFC <- lfcShrink(des, coef="genotype_kd_vs_wt", type = "apeglm")
  plotMA(resLFC)
  res_df <- resLFC %>% 
    as.data.frame() %>% 
    rownames_to_column("gene_symbol") %>% 
    as_tibble() %>% 
    mutate(qvalue = qvalue::qvalue(pvalue, fdr.level = 0.05)$qvalues)
  
  res_sig <- res_df %>% 
    filter(padj <= 0.05) %>% 
    arrange(desc(abs(log2FoldChange)))
  write_csv(res_sig, "data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_Results_signif.csv")

  saveRDS(dds, "data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_dds.rds")
  saveRDS(res_df, "data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_res_shrunk.rds")
}else{dds <- read_rds("data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_dds.rds")
      res_df <- read_rds("data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_res_shrunk.rds")}

# Volcano plot
res_vol <- res_df %>% 
  mutate(score = log2FoldChange*-log10(padj)) %>% 
  slice_max(score, n = 30) %>% 
  bind_rows(res_df %>% 
              mutate(score = log2FoldChange*-log10(padj)) %>% 
              slice_min(score, n = 30)) %>% 
  pull(gene_symbol)

EnhancedVolcano(res_df,
                lab = res_df$gene_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = paste0('Differential expression in aggregated agc1 samples vs aggregated control'),
                subtitle = "padj < 0.05 - |log2FC| thresholds at 1",
                pCutoff = 0.05, #0.05 cutoff
                FCcutoff = 1,
                axisLabSize = 14,
                titleLabSize = 14,
                subtitleLabSize = 10,
                captionLabSize = 10,
                pointSize = 2.5,
                labSize = 3.5,
                legendLabSize = 10,
                drawConnectors = T,
                caption = paste0("Significantly downregulated genes = ", 
                                 res_df %>% 
                                   filter(log2FoldChange < -1 & padj < 0.05) %>%
                                   nrow(),
                                 "     ",
                                 "Significantly upregulated genes = ",
                                 res_df %>% 
                                   filter(log2FoldChange > 1 & padj < 0.05) %>% 
                                   nrow()),
                selectLab = c("SREBF1", "SREBF2", "SLC25A13", "SLC25A12", res_vol)
) +
  theme(plot.caption = element_text(hjust = 0.5))

ggsave("plots/006a_RNASeq_patients_noTESSA/006a_4_Volcano_plot.png",
       h = 3500,
       w = 6000,
       dpi = 600,
       units = "px")

# Pathway analysis ----
fname <- "data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_mdf_hs.rds"
if(!file.exists(fname)){
  library(msigdbr) 
  msigdbr_species()
  msigdbr_collections() %>% print(n = Inf)
  mdf <- msigdbr(species = "Homo sapiens") %>%  
    filter(gs_cat %in% c("C2", "C5", "C7", "H"), 
           gs_subcat %in% c("GO:BP", "CP:WIKPATHWAYS", ""))
  saveRDS(mdf, file = fname)
}else(mdf <- readRDS(fname))

mlist <- split(x = mdf$gene_symbol, f = mdf$gs_name)

# Get signatures
fname <- "data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_fgsea_results_patients.rds"
if(!file.exists(fname)){
  
  library(fgsea)
  library(data.table)
  
  df <- res_df %>% as.data.frame() %>% 
    filter(!is.na(padj) & !is.na(log2FoldChange)) %>% 
    mutate(padj = replace(padj, padj == 0, 2.225074e-308)) 
  
  sig <- setNames(df$log2FoldChange, df$gene_symbol)
  #Use fgseaMultilevel for better accuracy than fgseaSimple (https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/fgsea/inst/doc/fgseaMultilevel-tutorial.html)
  fgseaRes <- fgseaMultilevel(pathways = mlist,
                              stats = sig,
                              eps = 0,
                              nPermSimple = 10000,
                              nproc = parallel::detectCores()-2,
  )  
  
  res_csv <- fgseaRes %>% 
    as_tibble() %>% 
    arrange(padj) %>% 
    rowwise() %>% 
    mutate(leadingEdge  = toString(leadingEdge)) %>% 
    ungroup()
  
  saveRDS(fgseaRes, fname)
  write_csv(res_csv, "data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_fgsea_results.csv")
  
}else{fgseaRes <- readRDS(fname)}

# Plots
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                      mlist, sig)

mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

png(paste0("plots/006a_RNASeq_patients_noTESSA/006a_07_agc1_collapsed_pathways.png"), h = 5000, w = 6500, res = 600)
plotGseaTable(mlist[mainPathways][1:20], sig, fgseaRes, 
              gseaParam = 0.5, colwidths = c(12,4,1.8,2,2))
dev.off()

## Barplot gsea
source("../../functions/p_star.R")
source("../../functions/pretty_path_label.R")
source("../../functions/plotEnrichment2.R")

path_int <- c("GOBP_NEUROPEPTIDE_SIGNALING_PATHWAY",
              "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_PLURIPOTENT_STEM_CELLS",
              "GOBP_LONG_TERM_SYNAPTIC_DEPRESSION",
              "GOBP_MONOAMINE_TRANSPORT",
              "REACTOME_GAP_JUNCTION_ASSEMBLY",
              "GOBP_DENDRITE_EXTENSION",
              "GOBP_POSITIVE_REGULATION_OF_DENDRITE_DEVELOPMENT",
              "GOBP_CELL_DIFFERENTIATION_IN_SPINAL_CORD",
              "REACTOME_NEURONAL_SYSTEM",
              "GOBP_POSITIVE_REGULATION_OF_LONG_TERM_SYNAPTIC_DEPRESSION",
              "GOBP_NEURON_DIFFERENTIATION",
              "GOBP_CENTRAL_NERVOUS_SYSTEM_DEVELOPMENT",
              "GOBP_FOREBRAIN_NEUROBLAST_DIVISION",
              "GOBP_EPIDERMAL_CELL_DIFFERENTIATION",
              "GOBP_METENCEPHALON_DEVELOPMENT",
              
              "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
              "REACTOME_MET_ACTIVATES_PTK2_SIGNALING",
              "GOBP_CONNECTIVE_TISSUE_DEVELOPMENT",
              "GOBP_SEMAPHORIN_PLEXIN_SIGNALING_PATHWAY_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE",
              "GOBP_DORSAL_ROOT_GANGLION_DEVELOPMENT",
              "GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
              "GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA",
              "GOBP_REGULATION_OF_EPITHELIAL_CELL_APOPTOTIC_PROCESS",
              "GOBP_EMBRYONIC_SKELETAL_SYSTEM_DEVELOPMENT",
              "REACTOME_SIGNALING_BY_PDGF",
              "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_RUNX2",
              "GOBP_MUSCLE_STRUCTURE_DEVELOPMENT",
              "GOBP_RESPONSE_TO_GROWTH_FACTOR",
              "GOBP_POSITIVE_REGULATION_OF_CELL_DEATH",
              "GOBP_MESENCHYMAL_CELL_DIFFERENTIATION")

ends <- fgseaRes %>% 
  filter(pathway %in% mainPathways) |> 
  na.omit() |> 
  # filter(pathway %in% path_int) %>%
  group_by(sign(NES)) |> 
  slice_max(abs(NES), n = 10) |> 
  ungroup() |> 
  arrange(NES) |> 
  mutate(across(pathway, ~str_replace(.x, "PLATELET_DERIVED_GROWTH_FACTOR_RECEPTOR", "PDGFR")))

png(paste0("plots/006a_RNASeq_patients_noTESSA/006a_08_agc1_patients_GSEA_barplot.png"), h = 4000, w = 5000, res = 600)

plt <- ggplot(ends, aes(x = c(1:nrow(ends)), y = NES, fill = as.factor(sign(-NES))), color = as.factor(sign(-NES))) +
  geom_bar(stat='identity', alpha = 0.85) +
  theme_light() +
  theme(legend.position = "none") +
  labs(title = "Best scoring pathways in AGC1 vs Control", y = "Combined Score", x="") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  coord_flip() +
  geom_text(aes(label = pretty_path_label(pathway) |> str_wrap(40), 
                y = ifelse(NES < 0, 0.05, -0.05),
                hjust = ifelse(NES < 0, 0, 1)),
            position = position_dodge(width = 0),
            size = 3.5,
            lineheight = 0.85) +
  geom_text(aes(label = p_star(padj)), 
            hjust = ifelse(sign(ends$NES) < 0, 1.3, -0.3),
            vjust = 0.75)


print(plt)
dev.off()

## Plot single GSEAs
pdf("plots/006a_RNASeq_patients_noTESSA/006a_09_agc1_patients_main_pathways_enrichment.pdf", w = 10, h = 5)
for(i in 1:length(path_int)){
  
  NES <- signif(fgseaRes[pathway == path_int[i]]$NES, 3)
  padj <- signif(fgseaRes[pathway == path_int[i]]$padj, 3)
  title <- paste0("GSEA in \n", 
                  pretty_path_label(path_int[i]) %>% str_wrap(50))
  
  p1 <- plotEnrichment2(mlist[[path_int[i]]], stats = sig) +
    labs(title = title, 
         subtitle = paste0("NES = ", NES, "  p.adj = ", padj)) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          plot.subtitle = element_text(size = 8, hjust = 0.5)) +
    xlim(0, 17000)
  
  ledge <- ends %>% 
    filter(pathway == path_int[i]) %>% 
    pull(leadingEdge) %>% 
    unlist() %>% 
    head(15)
  
  ledge_plot <- res_df %>%
    filter(gene_symbol %in% ledge) %>% 
    mutate(gene_symbol = fct_reorder(gene_symbol, log2FoldChange),
           fill_color = if_else(log2FoldChange < 0, "Neg", "Pos"))
  
  p2 <- ledge_plot %>% 
    ggplot(aes(gene_symbol, log2FoldChange, 
               fill = fill_color,
               color = fill_color)) +
    geom_bar(stat='identity', alpha = 0.7) +
    theme_light() +
    theme(legend.position = "none") +
    labs(title = "Leading edge genes expression", 
         y = "log2FoldChange", 
         x = "") +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 10)) +
    geom_text(aes(label = gene_symbol, 
                  y = ifelse(log2FoldChange < 0, 0.05, -0.05),
                  hjust = ifelse(log2FoldChange < 0, 1.1, -0.1)),
              position = position_dodge(width = 0),
              size = 3,
              lineheight = 0.85,
              color = "black") +
    geom_text(aes(label = p_star(padj), 
                  hjust = ifelse(sign(log2FoldChange) < 0, -0.1, 1.1),
                  vjust = 0.75),
              color = "black") +
    coord_flip() +
    scale_fill_manual(values = c("Pos" = "#F8766D", "Neg" = "#00BFC4")) +
    scale_color_manual(values = c("Pos" = "#F8766D", "Neg" = "#00BFC4")) 
  
  pf <- p1 + p2 +
    plot_layout(width = c(1.5,1))
  print(pf)
}
dev.off()

### Chromosome expression ----
# library(Organism.dplyr)
# library(org.Hs.eg.db)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# 
# # Add chromosome bands to deseq results
# supportedOrganisms() %>% print(n = Inf)
# src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
# 
# # keytypes(src)
# # map(keytypes(src), function(x){
# #   keys(src, keytype = x) %>% head()
# # })
# 
# # Create chromosomic map
# hs_chr <- select(src, 
#                  keys = res_df$gene_symbol,
#                  columns = c("alias", "cds_chrom", "gene_start"),
#                  keytype = "alias"
# ) %>% 
#   group_by(alias) %>% 
#   slice_head() %>% 
#   ungroup()
# 
# res_chr <- res_df %>% 
#   as.data.frame() %>% 
#   left_join(hs_chr %>% dplyr::select(alias, cds_chrom, gene_start), 
#             by = c("gene_symbol" = "alias")) %>% 
#   arrange(padj) %>% 
#   dplyr::select(gene_symbol, cds_chrom, gene_start, log2FoldChange, padj, pvalue) 
# 
# saveRDS(res_chr, "data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_DESeq_results_chromosome.rds")
# write_csv(res_chr , "data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_DESeq_results_chromosome.csv")
# 
# # Barplot chr coord
# library(ggrepel)
# 
# for(chr in str_subset(unique(res_chr$cds_chrom), "^chr")){
#   
#   subchr <- res_chr %>% 
#     filter(cds_chrom == chr) %>% 
#     filter(!is.na(log2FoldChange) & padj < 0.05) %>% 
#     mutate(start = factor(gene_start)) %>% 
#     rename("log2FC" = "log2FoldChange" )
#   
#   top_labels <- subchr %>% 
#     arrange(desc(abs(log2FC))) %>% 
#     slice_head(n = 20) %>% 
#     pull(gene_symbol)
#   
#   if(nrow(subchr) > 10){
#     
#     png(paste0("plots/006a_RNASeq_patients_noTESSA/006a_chromosome_expression/006a_11_chr", chr, "_expression.png"), h = 3000, w = 5000, res = 600)
#     p <- ggplot(subchr, aes(x = start, y = log2FC)) +
#       geom_bar(aes(fill = log2FC), stat = "identity") +
#       geom_smooth(aes(x = as.numeric(start)),
#                   #fill = "white",
#                   method = "loess", 
#                   formula = y ~ x,
#                   span = 0.4,
#                   lwd = 1.3) + 
#       labs(title = paste0(chr, " expression in aggregated ctr vs aggregated siSLC25A12 cells"),
#            subtitle = "padj < 0.05",) + 
#       xlab("Chromosomal coordinate") +
#       geom_label_repel(aes(label = ifelse(gene_symbol %in% top_labels, gene_symbol, "")),
#                        fontface = "bold",
#                        size = 3,
#                        color = "grey30",
#                        # segment.size = 1,
#                        # label.size = 1,
#                        box.padding = 0.5,
#                        max.overlaps = Inf,
#                        fill = "#00000020",
#                        seed = 41) +
#       scale_fill_viridis_c(aesthetics = c("fill", "segment.color", "color"), end = 0.9) +
#       scale_x_discrete(breaks = br <- levels(subchr$start)[floor(seq(1, nlevels(subchr$start), length.out = 10))],
#                        labels = scales::scientific(as.numeric(br))) +
#       theme(plot.title = element_text(hjust = 0.5),
#             plot.subtitle = element_text(hjust = 0.5)) +
#       theme_bw()
#     print(p)
#     dev.off()
#     
#   }
# }

# Dotplot log2fc ----
# Load automatically generated results and substitute the pre computed log2fc 
# with the one calculated not includin TESSA samples
auto_res_prim <- read_table("data/RNASeq_Lasorsa/BGI_F21FTSEUHT2056_HOMypeiR/Differentially_expressed_gene/Diff_exp/gene_diff.xls")
auto_res_orig <- auto_res_prim %>% 
  inner_join(res_df %>% select(gene_symbol, log2FoldChange, pvalue), by = "gene_symbol") %>% 
  mutate(`diffexp_log2fc_ctr-vs-agc` = log2FoldChange,
         `diffexp_deseq2_pvalue_ctr-vs-agc` = pvalue) %>% 
  select(-log2FoldChange, -pvalue)

auto_res <- auto_res_orig %>% 
  select(gene_symbol, contains("ctr-vs-agc")) %>% 
  rename_with(~gsub("diffexp_|deseq2_|_ctr-vs-agc", "", .x)) %>% 
  mutate(padj = p.adjust(pvalue, "BH"))

genes_fam <- c("^SLC2A", "^SLC38", "^LDH", "^GLS", "^MPC", "^SLC25", "^PDH", "^PDK",
               "^PDP", "^MDH", "^IDH", "^GDH", "GPT1|GPT2", "GOT1|GOT2", "^BCAT", 
               "ME1|ME2|ME3","ACLY", "PPARGC1A", "PPARGC1B", "^TFAM", "^CTH", 
               "^CDH", "^CAM", "^GPX", "^SOD", "^GRX", "^TXN", "^GSR", "^CAT", "^XDH", 
               "^ACOX", "^BCL", "^CASP", "^PARP", "^ANXA")

genes_plot <- map_dfr(genes_fam, ~filter(auto_res, str_detect(gene_symbol, .x))) %>% 
  filter(padj < 0.05, !is.na(padj))

targets <- genes_plot %>% 
  pull(gene_symbol) %>% 
  sort()

auto_new_l2fc <- auto_res_orig %>% 
  filter(gene_symbol %in% targets) %>%
  select(gene_symbol, contains("log2fc")) %>% 
  select(!contains(c("KO52", "TESSA"))) %>% 
  rename_with(~gsub("diffexp_log2fc_|diffexp_deseq2_", "", .x)) %>% 
  pivot_longer(!gene_symbol, names_to = "contrast", values_to = "log2FC")

auto_new_padj <- auto_res_orig %>% 
  filter(gene_symbol %in% targets) %>%
  select(gene_symbol, contains("pvalue")) %>% 
  select(!contains(c("KO52", "TESSA"))) %>% 
  rename_with(~gsub("diffexp_log2fc_|diffexp_deseq2_", "", .x)) %>% 
  pivot_longer(!gene_symbol, names_to = "contrast", values_to = "pvalue") %>% 
  group_by(contrast) %>% 
  mutate(padj = p.adjust(pvalue, "BH"),
         padj = ifelse(padj == 0, 1e-300, padj)) %>%
  ungroup()

all(rownames(auto_new_l2fc) == rownames(auto_new_padj)) # TRUE
all.equal(auto_new_l2fc$gene_symbol, auto_new_padj$gene_symbol) # TRUE

auto_new <- cbind(auto_new_l2fc, auto_new_padj %>% select(padj)) %>% 
  mutate(contrast = .$contrast %>% fct_relabel(~gsub("(MS|ATCC|ctr)(-vs-)(agc|agc_1|agc_2A|agc_2B)", "\\3\\2\\1", .x)),
         contrast = as_factor(contrast) %>% 
           fct_relevel(., "agc-vs-ATCC", after = Inf) %>% 
           fct_relevel(., "agc-vs-MS", after = Inf) %>%
           fct_relevel(., "agc-vs-ctr", after = Inf) %>% 
           fct_relevel(., rev),
         )

# for(i in c(1, seq(51, length(targets), 50))){
# 
#   auto_target <- targets[seq(i, i+49, 1)]
#   auto_target <- auto_target[!is.na(auto_target)]
#   
#   auto_new %>% 
#     filter(gene_symbol %in% auto_target, !is.na(padj)) %>% 
#     arrange(gene_symbol) %>% 
#     mutate(padj = ifelse(padj < 1e-100, 1e-100, padj)) %>% # Doing this because otherwise can't change the range of circles sizes using scale_size_continuous and large pvalues make all other supersmall
#     ggplot(aes(gene_symbol, contrast, size = -log10(padj))) +
#     geom_point(aes(fill = log2FC), shape = 21, color = "black") +
#     theme(axis.text.x = element_text(angle = 45)) +
#     labs(title = "Non aggregated and aggregated contrasts log2FoldChange",
#          subtitle = "Maximum padj = 0.05",
#          x = "Gene",
#          y = "Contrast") +
#     scale_color_distiller(palette = "Spectral", limits = c(-3,3), oob = scales::squish, aesthetics = "fill") +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     scale_size_continuous(labels = c(as.character(seq(0, 75, 25)), ">= 100"), breaks = c(0, 25, 50, 75, 100))
#   
#   ggsave(paste0("plots/006a_RNASeq_patients_noTESSA/dotplots_expression/006a_12_log2FC_heatmap_", i, "_", i+49, ".png"), 
#          h = 3000, w = 6200, dpi = 500, units = "px")
# }

# Oli Neu comparison
# res_mm <- readRDS("data/002_res.rds") %>% 
#   as.data.frame() %>% 
#   rownames_to_column("ensembl_gene") 
# 
# msdb <- msigdbr::msigdbr(species = "Mus musculus")
# 
# conv_table <- msdb %>% 
#   group_by(ensembl_gene) %>% 
#   arrange(desc(num_ortholog_sources)) %>% 
#   slice_head() %>% 
#   select(gene_symbol, entrez_gene, ensembl_gene, starts_with("human")) %>% 
#   ungroup()
# 
# res_join <- res_mm %>% 
#   select(ensembl_gene, log2FoldChange, padj, symbol) %>% 
#   as_tibble() %>% 
#   rename("gene_symbol" = "symbol") %>% 
#   left_join(conv_table %>% select(gene_symbol, human_gene_symbol), by = "gene_symbol") %>% 
#   inner_join(res_df %>% 
#                select(gene_symbol, log2FoldChange, padj) %>%
#                rename("human_gene_symbol" = gene_symbol), 
#              by = "human_gene_symbol", suffix = c("_mouse","_human"))
# 
# res_join %>% 
#   # filter(padj_mouse < 0.05 & padj_human < 0.05) %>% 
#   ggplot(aes(log2FoldChange_human, log2FoldChange_mouse)) +
#   ggpointdensity::geom_pointdensity(alpha = 0.7) + 
#   scale_color_viridis_c() +
#   geom_smooth(method = "lm", formula = "y ~ x") +
#   labs(title = "Correlation between human-derived siSLC25A12 cells and mouse OliNeu siSlc25a12 samples")
# ggsave("plots/006a_RNASeq_patients_noTESSA/006a_013_Correlation_with_OliNeu.png",
#        h = 1500,
#        w = 3000,
#        units = "px")

oli_sig <- readRDS("../Agc1_ATAC/data/007_OliNeu_signature.rds")
oli_ort <- orthologs(oli_sig$gene_symbol, species = "Mus musculus", human = F)
oli_sig_h <- inner_join(oli_sig, 
                        oli_ort |> 
                          select(human_symbol, symbol),
                        by = c("gene_symbol" = "symbol"))
merge_df <- oli_sig_h |> 
  inner_join(res_df |> 
               filter(padj < 0.05), 
             by = c("human_symbol" = "gene_symbol"))
