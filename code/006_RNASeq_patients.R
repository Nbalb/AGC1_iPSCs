library(EnhancedVolcano)
library(RColorBrewer)
library(ComplexHeatmap)
library(DESeq2)
library(msigdbr)
library(tidyverse)

# DESeq2 Differential Expression ----
if(!all(file.exists("data/RNASeq_Lasorsa/006_RNASeq_patients/006_res_shrunk.rds", 
                    "data/RNASeq_Lasorsa/006_RNASeq_patients/006_dds.rds"))){
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
                    cell_line = str_remove(sample, "_[1,2,3]$"))%>% 
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
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  png("plots/006_RNASeq_patients/006_1_Euclidean_distance.png", h = 3000, w = 4200, res = 600)
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
 
  ggsave("plots/006_RNAseq_patients/006_2_PCA_plot.png",
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
  
  png("plots/006_RNASeq_patients/006_1_Euclidean_distance_mean.png", h = 1200, w = 2000, res = 400)
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
  
  # Differential expression
  des <- DESeq(dds)
  res <- results(des)
  resLFC <- lfcShrink(dds, coef="genotype_kd_vs_wt", type = "apeglm")
  plotMA(resLFC)
  res_df <- resLFC %>% 
    as.data.frame() %>% 
    rownames_to_column("gene_symbol") %>% 
    as_tibble() %>% 
    mutate(qvalue = qvalue::qvalue(pvalue, fdr.level = 0.05)$qvalues)
  
  saveRDS(dds, "data/RNASeq_Lasorsa/006_RNASeq_patients/006_dds.rds")
  saveRDS(res_df, "data/RNASeq_Lasorsa/006_RNASeq_patients/006_res_shrunk.rds")
  }else{dds <- read_rds("data/RNASeq_Lasorsa/006_RNASeq_patients/006_dds.rds")
        res_df <- read_rds("data/RNASeq_Lasorsa/006_RNASeq_patients/006_res_shrunk.rds")}

# Load automatically generated results and check with our analysis
auto_res_orig <- read_table("data/RNASeq_Lasorsa/BGI_F21FTSEUHT2056_HOMypeiR/Differentially_expressed_gene/Diff_exp/gene_diff.xls")
auto_res <- auto_res_orig %>% 
  select(gene_symbol, contains("ctr-vs-agc")) %>% 
  rename_with(~gsub("diffexp_|deseq2_|_ctr-vs-agc", "", .x)) %>% 
  mutate(padj = p.adjust(pvalue, "BH"))

combo <- res_df %>%  
  select(-baseMean, -lfcSE) %>% 
  rename("log2fc" = "log2FoldChange") %>% 
  inner_join(auto_res, suffix = c("_our", "_auto"), by = "gene_symbol") 

combo %>% 
  filter(padj_our < 0.05, padj_auto < 0.05) %>% 
  ggplot(aes(log2fc_our, log2fc_auto)) +
  ggpointdensity::geom_pointdensity(size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  labs(title = "log2FC comparison between our analysis and automatic analysis",
       subtitle = "padj < 0.05")

ggsave("plots/006_RNAseq_patients/006_3_log2fc_comparison_our_vs_auto_padj.png",
       h = 3500,
       w = 5000,
       dpi = 600,
       units = "px")

# Volcano plot
EnhancedVolcano(res_df,
                lab = res_df$gene_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = paste0('Differential expression in aggregated agc1 samples vs aggregated control'),
                subtitle = "padj < 0.05",
                pCutoff = 0.05, #0.05 cutoff
                FCcutoff = 2,
                axisLabSize = 14,
                titleLabSize = 14,
                subtitleLabSize = 10,
                captionLabSize = 10,
                min.segment.length = 3,
                pointSize = 2.5,
                labSize = 3.5,
                legendLabSize = 10,
                drawConnectors = T,
                caption = paste0("Significantly downregulated genes = ", 
                                 res_df %>% 
                                   filter(log2FoldChange < -2 & padj < 0.05) %>%
                                   nrow(),
                                 "     ",
                                 "Significantly upregulated genes = ",
                                 res_df %>% 
                                   filter(log2FoldChange > 2 & padj < 0.05) %>% 
                                   nrow()),
                selectLab = c("SREBF1", "SREBF2", "SLC25A13", "SLC25A12",
                              res_df %>% 
                                filter(abs(log2FoldChange) > 5 | padj < 1e-150) %>% 
                                pull(gene_symbol))
) +
  theme(plot.caption = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(-10, -5, -2, 0, 2, 5, 10))

ggsave("plots/006_RNAseq_patients/006_4_Volcano_plot.png",
       h = 4000,
       w = 6500,
       dpi = 600,
       units = "px")

# Genes of interest ----
# Load tpm
auto_tpm_orig <- read_table("data/RNASeq_Lasorsa/BGI_F21FTSEUHT2056_HOMypeiR/Expression/gene_expression.xls")
auto_tpm <- auto_tpm_orig %>% 
  select(gene_symbol, starts_with("tpm"), -contains("KO"))

targets <-  c("Slc25a12", "Slc25a13", "Nat8l", "Aspa", "Acss1", "Srebf1",
              "Srebf2", "Fasn", "Olig1", "Olig2", "Olig3", "Sox10",
              "Pdgfra", "Cspg4", "Cldn11", "Mog", "Hdac9", "Hdac5", "Hdac2",
              "Hdac7", "Hdac3", "Hdac4", "Hdac1", "Hdac6", "Hdac11", "Hdac10", 
              "Hdac8", "Rai1" ,"Bcl2", "Bcl2l1", "Mycn", "Psmb8") %>% 
  toupper()

all(targets %in% res_df$gene_symbol)

auto_tpm_long <- auto_tpm %>% 
  filter(gene_symbol %in% targets) %>%
  left_join(auto_res %>% select(gene_symbol, log2fc, padj), by = "gene_symbol") %>% 
  arrange(gene_symbol) %>%
  # Need to set the labels for the facets, and in order to change the labels and not the level of the factor
  # its best to use fct_relabel. The best option is probably to just change gene_symbol and use it as the 
  # facetting variable. But I'll just leave it like this in case I ever need an example of fct_relabel()
    mutate(gene_label = paste0(gene_symbol, "log2FC = ", .$log2fc, "padj = ", .$padj),
           gene_symbol = as_factor(gene_symbol),
           gene_symbol = fct_relabel(gene_symbol, function(x) paste0(gene_symbol, 
                                                                     "\nlog2FC = ", .$log2fc %>% scales::scientific(), 
                                                                     " - padj = ", .$padj %>% scales::scientific()))) %>% 
  pivot_longer(starts_with("tpm"), names_to = "cell_line") %>% 
  mutate(cell_label = str_extract(cell_line, "agc_[0-9A-Z]{1,2}|ctr[A-Z]{2,4}") %>% as_factor(),
         bg_color = ifelse(log2fc < 0 & padj < 0.05, "lightblue",
                           ifelse(log2fc > 0 & padj < 0.05, "red", "grey30")))

auto_tpm_long %>% 
  slice_head(prop = 0.5) %>% 
ggplot(aes(cell_label, value, color = cell_label)) +
  geom_tile(aes(xmin = -Inf, x = Inf, ymin = -Inf, y = Inf, fill = bg_color), alpha = 0.2, inherit.aes = F) +
  geom_point(size = 3, alpha = 0.7) +
  facet_wrap(~ gene_symbol, nrow = 4, scales = "free_y") + 
  scale_color_viridis_d(end = 0.8) +
  scale_y_continuous(expand = c(.15, .15)) +
  labs(title = "TPM counts for relevant genes",
       x = "Cell Line",
       y = "TPM",
       fill = "Cell Line") +
  scale_fill_manual("Important aspect", values=c("lightblue","grey","red"))
  theme(axis.text.x = element_text(size = 7),
        strip.text.x = element_text(size = 8))

ggsave("plots/006_RNAseq_patients/006_5_Relevant_genes_tpm1.png",
       h = 4320,
       w = 7680,
       units = "px",
       dpi = 600)

auto_tpm_long %>% 
  slice_tail(prop = 0.5) %>% 
  ggplot(aes(cell_label, value, color = cell_label)) +
  geom_point(size = 3, alpha = 0.7) +
  facet_wrap(~ gene_symbol, nrow = 4, scales = "free_y") + 
  scale_color_viridis_d(end = 0.8) +
  scale_y_continuous(expand = c(.15, .15)) +
  labs(title = "TPM counts for relevant genes",
       x = "Cell Line",
       y = "TPM",
       fill = "Cell Line") +
  theme(axis.text.x = element_text(size = 7),
        strip.text.x = element_text(size = 8))

ggsave("plots/006_RNAseq_patients/006_6_Relevant_genes_tpm2.png",
       h = 4320,
       w = 7680,
       units = "px",
       dpi = 600)

auto_tpm %>% 
  pivot_longer(cols = starts_with("tpm"), names_to = "sample", values_to = "TPM") %>% 
  mutate(sample = str_remove(sample, "tpm_")) %>% 
  ggplot(aes(gene_symbol, sample, size = TPM)) +  
  geom_point(alpha = 0.7, color = "midnightblue") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(title = "Relevant genes expression")


# Pathway analysis ----
fname <- "data/RNASeq_Lasorsa/006_RNASeq_patients/006_mdf_hs.rds"
if(!file.exists(fname)){
  library(msigdbr) 
  msigdbr_species()
  msigdbr_collections() %>% print(n = Inf)
  mdf <- msigdbr(species = "Homo sapiens") %>%     # Retrieve all Dm gene sets
    filter(gs_cat %in% c("C1", "C2", "C5", "C7", "H"), 
           gs_subcat %in% c("", "GO:MF", "GO:BP", "GO:CC", "CP:KEGG", "CP:WIKIPATHWAYS", "CP:REACTOME"))
  saveRDS(mdf, file = fname)
}else(mdf <- readRDS(fname))

mlist <- split(x = mdf$gene_symbol, f = mdf$gs_name)

# Get signatures
fname <- "data/RNASeq_Lasorsa/006_RNASeq_patients/006_fgsea_results_patients.rds"
if(!file.exists(fname)){
  
  library(fgsea)
  library(data.table)
  
  df <- auto_res %>% as.data.frame() %>% 
    filter(!is.na(padj) & !is.na(log2fc)) %>% 
    mutate(padj = replace(padj, padj == 0, 2.225074e-308)) 
  
  sig <- setNames(df$log2fc, df$gene_symbol)
  #Use fgseaMultilevel for better accuracy than fgseaSimple (https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/fgsea/inst/doc/fgseaMultilevel-tutorial.html)
  fgseaRes <- fgseaMultilevel(pathways = mlist,
                              stats = sig,
                              eps = 0,
                              nPermSimple = 10000,
                              nproc = parallel::detectCores()-2,
  )  
  
  saveRDS(fgseaRes, fname)
  write_csv(fgseaRes[padj < 0.05][order(padj)], "data/RNASeq_Lasorsa/006_RNASeq_patients/006_fgsea_results_significant.csv")
  
}else{fgseaRes <- readRDS(fname)}

# Plots
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], 
                                      mlist, sig)

mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

png(paste0("plots/006_RNASeq_patients/006_07_agc1_collapsed_pathways.png"), h = 5000, w = 6500, res = 600)
plotGseaTable(mlist[mainPathways][1:20], sig, fgseaRes, 
              gseaParam = 0.5, colwidths = c(12,4,1.8,2,2))
dev.off()

## Barplot gsea
source("../../functions/p_star.R")
source("../../functions/pretty_path_label.R")
source("../../functions/plotEnrichment2.R")

ends <- fgseaRes %>%
  filter(pathway %in% mainPathways) %>% 
  filter(NES > 0 & padj < 0.05) %>% 
  slice_min(padj, n = 10) %>% 
  arrange(-NES) %>% 
  bind_rows(fgseaRes %>%
              filter(NES < 0 & padj < 0.05) %>% 
              slice_min(padj, n = 10) %>% 
              arrange(-NES)) %>% 
  arrange(NES)


png(paste0("plots/006_RNAseq_patients/006_08_agc1_patients_GSEA_barplot.png"), h = 3500, w = 4500, res = 600)

plt <- ggplot(ends, aes(x = c(1:nrow(ends)), y = NES, fill = as.factor(sign(-NES)))) +
  geom_bar(stat='identity') +
  theme_light() +
  theme(legend.position = "none") +
  labs(title = "Best scoring pathways in siAgc1 OliNeu cells", y = "Combined Score", x="") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  coord_flip() +
  geom_text(aes(label = pretty_path_label(pathway), 
                y = ifelse(NES < 0, 0.05, -0.05),
                hjust = ifelse(NES < 0, 0, 1)),
            position = position_dodge(width = 0),
            size = 3,
            lineheight = 0.85) +
  geom_text(aes(label = p_star(padj)), 
            hjust = ifelse(sign(ends$NES) < 0, 1.3, -0.3),
            vjust = 0.75)


print(plt)
dev.off()

## Plot single GSEAs
pdf(paste0("plots/006_RNAseq_patients/006_09_agc1_patients_main_pathways_enrichment.pdf"), w = 5, h = 3)
for(o in 1:20){
  
  NES <- signif(fgseaRes[pathway == mainPathways[o]]$NES, 3)
  padj <- signif(fgseaRes[pathway == mainPathways[o]]$padj, 3)
  title <- paste0("GSEA in \n", 
                  pretty_path_label(mainPathways[o]) %>% str_wrap(50))
  
  p <- plotEnrichment2(mlist[[mainPathways[o]]], stats = sig) +
    labs(title = title, 
         subtitle = paste0("NES = ", NES, "  p.adj = ", padj)) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          plot.subtitle = element_text(size = 8, hjust = 0.5)) +
    xlim(0, 10500)
  print(p)
}
dev.off()

### Chromosome expression ----
library(Organism.dplyr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Add chromosome bands to deseq results
supportedOrganisms() %>% print(n = Inf)
src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")

keytypes(src)
map(keytypes(src), function(x){
  keys(src, keytype = x) %>% head()
})
# Create chromosomic map
hs_chr <- select(src, 
                 keys = res_df$gene_symbol,
                 columns = c("alias", "cds_chrom", "gene_start"),
                 keytype = "alias"
) %>% 
  group_by(alias) %>% 
  slice_head() %>% 
  ungroup()

res_chr <- res_df %>% 
  as.data.frame() %>% 
  left_join(hs_chr %>% dplyr::select(alias, cds_chrom, gene_start), 
            by = c("gene_symbol" = "alias")) %>% 
  arrange(padj) %>% 
  dplyr::select(gene_symbol, cds_chrom, gene_start, log2FoldChange, padj, pvalue) 

saveRDS(res_chr, "data/RNASeq_Lasorsa/006_RNASeq_patients/006_DESeq_results_chromosome.rds")
write_csv(res_chr, "data/RNASeq_Lasorsa/006_RNASeq_patients/006_DESeq_results_chromosome.csv")

# Barplot best DE genes
source("../../functions/p_star.R")

de_plot <- res_chr %>%
  mutate(gene_label = paste0(gene_symbol, " - ", cds_chrom)) %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  slice_max(n = 10, order_by = log2FoldChange) %>%
  arrange(desc(log2FoldChange)) %>% 
  bind_rows(res_chr %>% 
              mutate(gene_label = paste0(gene_symbol, " - ", cds_chrom)) %>% 
              filter(padj < 0.05 & log2FoldChange < 0) %>% 
              slice_min(n = 10, order_by = log2FoldChange) %>%
              arrange(desc(log2FoldChange))) %>% 
  mutate(gene_label = fct_reorder(gene_label, log2FoldChange)) 

de_plot %>% 
ggplot(aes(gene_label, log2FoldChange, 
           fill = as_factor(-1*sign(log2FoldChange)),
           color = as_factor(-1*sign(log2FoldChange)))) +
  geom_bar(stat='identity', alpha = 0.7) +
  theme_light() +
  theme(legend.position = "none") +
  labs(title = "Best DE genes in aggregated siSLC25A12 vs ctr", y = "log2FoldChange", x="") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10)) +
  coord_flip() +
  geom_text(aes(label = gene_label, 
                y = ifelse(log2FoldChange < 0, 0.05, -0.05),
                hjust = ifelse(log2FoldChange < 0, 0, 1)),
            position = position_dodge(width = 0),
            size = 3,
            lineheight = 0.85,
            color = "black") +
  geom_text(aes(label = p_star(padj), 
            hjust = ifelse(sign(log2FoldChange) < 0, 1.3, -0.3),
            vjust = 0.75),
            color = "black")
ggsave("plots/006_RNAseq_patients/006_10_Best_DE_genes_barplot.png",
        h = 3500,
        w = 4500,
        units = "px",
        dpi = 600)

# Barplot chr coord
library(ggrepel)

for(chr in str_subset(unique(res_chr$cds_chrom), "^chr")){
  
  subchr <- res_chr %>% 
    filter(cds_chrom == chr) %>% 
    filter(!is.na(log2FoldChange) & padj < 0.05) %>% 
    mutate(start = factor(gene_start)) %>% 
    rename("log2FC" = "log2FoldChange" )
  
  top_labels <- subchr %>% 
    arrange(desc(abs(log2FC))) %>% 
    slice_head(n = 20) %>% 
    pull(gene_symbol)
  
  if(nrow(subchr) > 10){
    
    png(paste0("plots/006_RNASeq_patients/006_chromosome_expression/006_11_chr", chr, "_expression.png"), h = 3000, w = 5000, res = 600)
    p <- ggplot(subchr, aes(x = start, y = log2FC)) +
      geom_bar(aes(fill = log2FC), stat = "identity") +
      geom_smooth(aes(x = as.numeric(start)),
                  #fill = "white",
                  method = "loess", 
                  formula = y ~ x,
                  span = 0.4,
                  lwd = 1.3) + 
      labs(title = paste0(chr, " expression in aggregated ctr vs aggregated siSLC25A12 cells"),
           subtitle = "padj < 0.05",) + 
      xlab("Chromosomal coordinate") +
      geom_label_repel(aes(label = ifelse(gene_symbol %in% top_labels, gene_symbol, ""),
                           color = log2FC,
                           segment.color = log2FC),
                       fontface = "bold",
                       size = 3,
                       # segment.size = 1,
                       # label.size = 1,
                       box.padding = 0.5,
                       max.overlaps = Inf,
                       fill = "#00000020",
                       seed = 41) +
      scale_fill_viridis_c(aesthetics = c("fill", "segment.color", "color"), end = 0.9) +
      scale_x_discrete(breaks = br <- levels(subchr$start)[floor(seq(1, nlevels(subchr$start), length.out = 10))],
                       labels = scales::scientific(as.numeric(br))) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)) +
      theme_bw()
    print(p)
    dev.off()
    
  }
}

# Heatmap log2fc
auto_new_l2fc <- auto_res_orig %>% 
  filter(gene_symbol %in% targets) %>%
  select(gene_symbol, contains("log2fc")) %>% 
  select(!contains(c("KO52"))) %>% 
  rename_with(~gsub("diffexp_log2fc_|diffexp_deseq2_", "", .x)) %>% 
  pivot_longer(!gene_symbol, names_to = "contrast", values_to = "log2FC")

auto_new_padj <- auto_res_orig %>% 
  filter(gene_symbol %in% targets) %>%
  select(gene_symbol, contains("pvalue")) %>% 
  select(!contains(c("KO52"))) %>% 
  rename_with(~gsub("diffexp_log2fc_|diffexp_deseq2_", "", .x)) %>% 
  pivot_longer(!gene_symbol, names_to = "contrast", values_to = "pvalue") %>% 
  group_by(contrast) %>% 
  mutate(padj = p.adjust(pvalue, "BH"),
         padj = ifelse(padj == 0, 1e-300, padj)) %>%
  ungroup()

all(rownames(auto_new_l2fc) == rownames(auto_new_padj)) # TRUE
all.equal(auto_new_l2fc$gene_symbol, auto_new_padj$gene_symbol) # TRUE

auto_new <- cbind(auto_new_l2fc, auto_new_padj %>% select(padj)) %>% 
  mutate(contrast = .$contrast %>% fct_relabel(~gsub("(TESSA|MS|ATCC|ctr)(-vs-)(agc|agc_1|agc_2A|agc_2B)", "\\3\\2\\1", .x)),
         contrast = as_factor(contrast) %>% 
           fct_relevel(., "agc-vs-ATCC", after = Inf) %>% 
           fct_relevel(., "agc-vs-MS", after = Inf) %>% 
           fct_relevel(., "agc-vs-TESSA", after = Inf) %>% 
           fct_relevel(., "agc-vs-ctr", after = Inf) %>% 
           fct_relevel(., rev))

auto_new %>% 
  mutate(padj = ifelse(padj < 1e-100, 1e-100, padj)) %>% # Doing this because otherwise can't change the range of circles sizes using scale_size_continuous and large pvalues make all other supersmall
  ggplot(aes(gene_symbol, contrast, size = -log10(padj))) +
  geom_point(aes(fill = log2FC), shape = 21, color = "black") +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(title = "Non aggregated and aggregated contrasts log2FoldChange") +
  scale_color_distiller(palette = "Spectral", limits = c(-3,3), oob = scales::squish, aesthetics = "fill") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(labels = c(as.character(seq(0, 75, 25)), ">= 100"))
  
ggsave("plots/006_RNASeq_patients/006_12_log2FC_heatmap.png", h = 4000, w = 6000, dpi = 600, units = "px")

# Oli Neu comparison
res_mm <- readRDS("data/002_res.rds") %>% 
  as.data.frame() %>% 
  rownames_to_column("ensembl_gene") 

msdb <- msigdbr::msigdbr(species = "Mus musculus")

conv_table <- msdb %>% 
  group_by(ensembl_gene) %>% 
  arrange(desc(num_ortholog_sources)) %>% 
  slice_head() %>% 
  select(gene_symbol, entrez_gene, ensembl_gene, starts_with("human")) %>% 
  ungroup()

res_join <- res_mm %>% 
  select(ensembl_gene, log2FoldChange, padj, symbol) %>% 
  rename("gene_symbol" = symbol) %>% 
  left_join(conv_table %>% select(gene_symbol, human_gene_symbol), by = "gene_symbol") %>% 
  inner_join(res_df %>% 
               select(gene_symbol, log2FoldChange, padj) %>%
               rename("human_gene_symbol" = gene_symbol), 
             by = "human_gene_symbol", suffix = c("_mouse","_human"))
  
res_join %>% 
  # filter(padj_mouse < 0.05 & padj_human < 0.05) %>% 
  ggplot(aes(log2FoldChange_human, log2FoldChange_mouse)) +
  ggpointdensity::geom_pointdensity(alpha = 0.7) + 
  scale_color_viridis_c() +
  geom_smooth(method = "lm", formula = "y ~ x") +
  labs(title = "Correlation between human-derived siSLC25A12 cells and mouse OliNeu siSlc25a12 samples")
ggsave("plots/006_RNAseq_patients/006_013_Correlation_with_OliNeu.png",
       h = 1500,
       w = 3000,
       units = "px")
