library(tidyverse)
library(fgsea)
library(data.table)

#Look for metabolic pathways
res_df <- read_rds("data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_res_shrunk.rds")

fname <- "data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_mdf_hs.rds"
if(!file.exists(fname)){
  library(msigdbr) 
  msigdbr_species()
  msigdbr_collections() %>% print(n = Inf)
  mdf <- msigdbr(species = "Homo sapiens") %>%  
    filter(gs_cat %in% c("C1", "C2", "C5", "C7", "H"), 
           gs_subcat %in% c("CP:WIKIPATHWAYS", "GO:BP"))
  saveRDS(mdf, file = fname)
}else(mdf <- readRDS(fname))

mdf %>% 
  filter(gs_subcat == "CP:REACTOME", 
         str_detect(gs_name, "HISTONE")) %>% 
  group_by(gs_name) %>% 
  mutate(n = n()) %>% 
  slice_head() %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  select(gs_name, n) %>% 
  print(n = Inf)

metab_entr <- c("GLUCOSE|GLYCOLYSIS|PENTOSE|GLUTAMATE|
                               PYRUVATE|MITOCHONDRIA|CITRIC|RESPIRATORY|
                               KETONE|FATTY_ACIDS|UREA|ADHESION|CATENIN|
                               OXIDATIVE_STRESS|NADPH|PURINE|
                               PYRIMIDINE|APOPTOSIS|CASPASE|HISTONE")

metab_list <- mdf %>% 
  filter(gs_subcat == "CP:WIKIPATHWAYS")

mlist <- split(x = metab_list$gene_symbol, f = metab_list$gs_name)

# Get signatures
fname <- "data/RNASeq_Lasorsa/006b/006b_fgsea_results_metabolic_pathways_patients.rds"
if(!file.exists(fname)){

  df <- res_df %>% as.data.frame() %>% 
    filter(!is.na(padj) & !is.na(log2FoldChange)) %>% 
    mutate(padj = replace(padj, padj == 0, 2.225074e-308)) 
  
  sig <- setNames(df$log2FoldChange, df$gene_symbol)
  #Use fgseaMultilevel for better accuracy than fgseaSimple (https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/fgsea/inst/doc/fgseaMultilevel-tutorial.html)
  fgseaRes <- fgseaMultilevel(pathways = mlist,
                               stats = sig,
                               eps = 0,
                               nPermSimple = 10000,
                               nproc = parallel::detectCores()-1,
  )  
  
  saveRDS(fgseaRes, fname)
  write_csv(fgseaRes[order(padj)], "data/RNASeq_Lasorsa/006b/006b_fgsea_results_significant.csv")
  
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

ends <- fgseaRes %>%
  filter(pathway %in% path_int) %>%
  arrange(NES)


png(paste0("plots/006a_RNASeq_patients_noTESSA/006a_08_agc1_patients_GSEA_barplot.png"), h = 3500, w = 4500, res = 600)

plt <- ggplot(ends, aes(x = c(1:nrow(ends)), y = NES, fill = as.factor(sign(-NES))), color = as.factor(sign(-NES))) +
  geom_bar(stat='identity', alpha = 0.85) +
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
pdf("plots/006a_RNASeq_patients_noTESSA/006a_09_agc1_patients_main_pathways_enrichment.pdf", w = 10, h = 5)
for(o in 1:length(path_int)){
  
  NES <- signif(fgseaRes[pathway == path_int[o]]$NES, 3)
  padj <- signif(fgseaRes[pathway == path_int[o]]$padj, 3)
  title <- paste0("GSEA in \n", 
                  pretty_path_label(path_int[o]) %>% str_wrap(50))
  
  p1 <- plotEnrichment2(mlist[[path_int[o]]], stats = sig) +
    labs(title = title, 
         subtitle = paste0("NES = ", NES, "  p.adj = ", padj)) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          plot.subtitle = element_text(size = 8, hjust = 0.5)) +
    xlim(0, 17000)
  
  ledge <- ends %>% 
    filter(pathway == path_int[o]) %>% 
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
