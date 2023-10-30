library(EnhancedVolcano)
library(RColorBrewer)
library(ComplexHeatmap)
library(msigdbr)
library(gggsea)
source("funcs/aREA_GSEA_Plot.R")
library(tidyverse)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")
conflicted::conflict_prefer("slice", "dplyr")

# Load pathways and results
pathways <- read_csv("results/pathways_to_search_Lasorsa.csv", col_names = c("name", "x2", "gene_symbol")) %>% 
  select(1, 3)
msig <- msigdbr()
hif_pathways <- msig %>% 
  filter(str_detect(gs_name, "HIF|KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM")) %>% 
  select(gs_name, gene_symbol, gs_description) %>% 
  rename(name = gs_name)

res_df <- read_rds("data/RNASeq_Lasorsa/006a_RNASeq_patients_noTESSA/006a_res_shrunk.rds")

# Load automatically generated results and substitute the pre computed log2fc 
# with the one calculated excluding TESSA samples
notessa_plots <- function(pathways, plot_dir = "plots/006d/", row_fontsize = 8){
    
  pathway_un <- unique(pathways$name)
  targets <- pathways %>% 
    pull(gene_symbol)
  dir.create(plot_dir) %>% suppressWarnings()
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
           padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
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
  
  if(!file.exists("data/006d/results_noTessa.rds")){
    saveRDS(auto_new, "data/006d/results_noTessa.rds")
  }
  
  # Heatmaps of top differentially expressed genes
  best_diff <- auto_new %>% 
    filter(contrast == "agc-vs-ctr",
           !is.na(padj)) %>% 
    mutate(stat = sign(log2FC)*-log10(padj)) %>% 
    arrange(-stat) %>% 
    slice(-(26:(n()-25))) 
  
  auto_tpm_orig <- read_table("data/RNASeq_Lasorsa/BGI_F21FTSEUHT2056_HOMypeiR/Expression/gene_expression.xls")
  auto_tpm <- auto_tpm_orig %>% 
    select(gene_symbol, starts_with("tpm"), -contains("KO"), -contains("TESSA")) 
  
  best_diff_tpm <- auto_tpm %>%
    filter(gene_symbol %in% best_diff$gene_symbol) %>% 
    arrange(match(gene_symbol, best_diff$gene_symbol)) %>% 
    mutate(across(where(is.numeric), ~ log10(.x + 1))) %>% 
    rename_with(~ gsub("tpm_", "", .x)) %>% 
    relocate(starts_with("ctr")) %>% 
    column_to_rownames("gene_symbol") %>% 
    as.matrix() %>% 
    t() %>% 
    scale() %>% 
    t()
  
  colors_hm <- circlize::colorRamp2(
    breaks = seq(-1, 1, length = 3), 
    colors = viridis::plasma(n = 3), 
    space = "LAB")
  
  png(paste0(plot_dir, "/001_Heatmap_best_genes.png"), h = 3000, w = 4200, res = 400)
  Heatmap(best_diff_tpm,
          cluster_columns = F,
          col = colors_hm,
          column_labels = str_replace_all(colnames(best_diff_tpm), "_", " "),
          row_names_gp = grid::gpar(fontsize = 8),
          column_names_rot = 30,
          column_split = c(rep(1, 6), rep(2, 9)),
          column_title = "Top and bottom best differentially expressed genes",
          heatmap_legend_param = list(
            title = "log10-scaled \nTPM")
  )
  dev.off()
  
  
  # Pathways heatmaps ----
  dir.create(paste0(plot_dir, "/Heatmaps_pathway/")) %>% suppressWarnings()
  for(i in pathway_un){
    
    path_genes <- pathways %>% 
      filter(name == i) %>% 
      pull(gene_symbol)
    
    path_mtx <- auto_tpm %>%
      filter(gene_symbol %in% path_genes) %>% 
      #arrange(as.numeric(str_extract(gene_symbol, "[0-9]{1,2}$"))) %>% 
      mutate(across(where(is.numeric), ~ log10(.x + 1))) %>% 
      filter(!if_all(where(is.numeric), ~ .x == 0)) %>% 
      rename_with(~ gsub("tpm_", "", .x)) %>% 
      relocate(starts_with("ctr")) %>% 
      column_to_rownames("gene_symbol") %>% 
      as.matrix() %>% 
      t() %>% 
      scale() %>% 
      t()
    
    # if(nrow(path_mtx) > 30){
    #   row_fontsize <- 8
    #   plot_h <- 3000
    #   plot_w <- 3000
    # }else{row_fontsize <- 8
    # plot_h <- 1500
    # plot_w <- 2500}
    plot_h <- 2000+(max(0, nrow(path_mtx)-30))*50
    plot_w <- 2500
    
    colors_hm <- circlize::colorRamp2(
      breaks = seq(-1, 1, length = 3), 
      colors = viridis::plasma(n = 3), 
      space = "LAB")
    
    pathway_file_name <- str_replace_all(i, "[:punct:]|\\s", "_")
    png(paste0(plot_dir, "/Heatmaps_pathway/", pathway_file_name, ".png"), h = plot_h, w = plot_w, res = 400)
    hm <- Heatmap(path_mtx,
                  cluster_columns = F,
                  col = colors_hm,
                  column_labels = str_replace_all(colnames(best_diff_tpm), "_", " "),
                  row_names_gp = grid::gpar(fontsize = row_fontsize),
                  column_names_gp = grid::gpar(fontsize = 8),
                  column_names_rot = 30,
                  column_split = c(rep(1, 6), rep(2, 9)),
                  column_title = paste0(str_replace_all(i, "_", " ") %>% str_to_title(), " genes expression"),
                  heatmap_legend_param = list(
                    title = "log10-scaled \nTPM")
    )
    print(hm)
    dev.off()
  }
  
  # Dotplot pathways ----
  dir.create(paste0(plot_dir, "/dotplots_expression/")) %>% suppressWarnings()
  for(i in pathway_un){
    
    pathway_sel <- pathways %>% 
      filter(name == i)
   
    auto_target <- pathway_sel %>% 
      pull(gene_symbol)
    
    dot_w <- 3600+(max(0, length(auto_target)-30))*100
    subtitle <- if(!is.null(pathways$gs_description) %>% suppressWarnings()){
      paste0(pathway_sel$gs_description[1], "\nOnly points with padj < 0.05 are represented")
      
    }else{"Only points with padj < 0.05 are represented"}
    
    auto_new %>% 
      filter(gene_symbol %in% auto_target, !is.na(padj)) %>% 
      mutate(gene_symbol = fct_reorder(gene_symbol, log2FC)) %>% 
      mutate(logpadj = -log10(padj),
             logpadj = ifelse(logpadj >= 20, 20, logpadj)) %>% # Doing this because otherwise can't change the range of circles sizes using scale_size_continuous and large pvalues make all other supersmall
      ggplot(aes(gene_symbol, contrast, size = logpadj)) +
      geom_point(aes(fill = log2FC), shape = 21, color = "black") +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7)) +
      labs(title = paste(i, " - log2FoldChange"),
           subtitle = subtitle,
           x = "Gene",
           y = "Contrast") +
      scale_color_distiller(palette = "Spectral", limits = c(-2,2), oob = scales::squish, aesthetics = "fill") +
      scale_size_continuous(limits = c(-log10(0.05), 20), # This guarantees that points with padj > 0.05 become a blank space
                            range = c(1, 5.5),
                            breaks = seq(0, 20, 5),
                            labels = c(as.character(seq(0, 15, 5)), "20+"),
                            name = expression(paste(-log[10],"padj"))
      )
    
    ggsave(paste0(plot_dir, "/dotplots_expression/006c_log2FC_dotplot_", str_replace_all(i, "[:punct:]|\\s", "_"), ".png"), 
           h = 2000, w = dot_w, dpi = 400, units = "px", limitsize = F)
  }
}

notessa_plots(pathways)
notessa_plots(hif_pathways, "plots/006d/HIF")

# New pathways GSEA ----
library(fgsea)

path_gsea <- hif_pathways %>% 
  mutate(name = str_remove_all(name, " - [0-9]"))

mlist <- split(path_gsea$gene_symbol, path_gsea$name)
stat <- res_df %>% 
  filter(!is.na(padj)) %>% 
  mutate(padj = ifelse(padj == 0, .Machine$double.xmin, padj),
         stat = sign(log2FoldChange)*-log10(padj)) %>% 
  drop_na(stat) %>% 
  arrange(desc(stat)) %>% 
  pull(stat, gene_symbol) 

set.seed(41)
fgsea_res <- fgseaMultilevel(pathways = mlist,
                             stats = stat,
                             nproc = 6)

# source("../../functions/plotEnrichment2.R")
# png("plots/006d/GSEA_plot_ATP_synthase.png", h = 1500, w = 3000, res = 400)
# plotEnrichment2(mlist[["Respiratory chain/ATP synthase"]], stat)+
#   labs(title = "Respiratory chain/ATP synthase", 
#        subtitle = paste0("NES = ", round(NES, 3), "  p.adj = ", scales::scientific(padj))) +
#   theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
#         plot.subtitle = element_text(size = 8, hjust = 0.5)) +
#   xlim(0,20000)
# dev.off()

# Check HIF pathways
hif_list <- split(hif_pathways$gene_symbol, hif_pathways$name)
set.seed(41)
hif_gsea_res <- fgseaMultilevel(pathways = hif_list,
                                stats = stat)
gsea_signif <- hif_gsea_res[padj < 0.05]
hif_list_plot <- hif_list[gsea_signif$pathway]

df <- gseaCurve(stat, hif_list_plot, gsea_signif)
ggplot2::ggplot() + 
  geom_gsea(df) + 
  labs(title = "HIF-related pathways GSEA") +
  theme_gsea(7) +
  theme(plot.title = element_text(size = 18, face = "bold"))
ggsave("plots/006d/HIF/gsea_signif.png", h = 2000, w = 3000, units = "px")  

# aREA Plot
elv_merge <- c(setNames(rep(-1, length(hif_list_plot[["ELVIDGE_HIF1A_AND_HIF2A_TARGETS_DN"]])),
                        hif_list_plot[["ELVIDGE_HIF1A_AND_HIF2A_TARGETS_DN"]]),
               setNames(rep(1, length(hif_list_plot[["ELVIDGE_HIF1A_AND_HIF2A_TARGETS_UP"]])),
                        hif_list_plot[["ELVIDGE_HIF1A_AND_HIF2A_TARGETS_UP"]]))
aplot <- aREA_GSEA_Plot(current.geneset = elv_merge, 
               current.sig = stat,
               main.title = "HIF1A and HIF2A targets GSEA")
ggsave(plot = aplot, "plots/006d/HIF/area_plot_elvidge_hif1a_hif2a.png", h = 2000, w = 3000, units = "px")

hif1_merge <- c(setNames(rep(-1, length(hif_list[["ELVIDGE_HIF1A_TARGETS_DN"]])),
                          hif_list[["ELVIDGE_HIF1A_TARGETS_DN"]]),
                 setNames(rep(1, length(hif_list[["ELVIDGE_HIF1A_TARGETS_UP"]])),
                          hif_list[["ELVIDGE_HIF1A_TARGETS_UP"]]))
hif_plot <- aREA_GSEA_Plot(current.geneset = hif1_merge,
                           current.sig = stat,
                           main.title = "HIF1A targets GSEA")
ggsave(plot = hif_plot, "plots/006d/HIF/area_plot_elvidge_hif1a.png", h = 2000, w = 3000, units = "px")
