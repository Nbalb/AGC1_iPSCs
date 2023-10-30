library(ComplexHeatmap)
library(tidyverse)
conflicted::conflict_prefer_all("dplyr", quiet = T)

# Heatmaps of top differentially expressed
res_df <- read_rds("data/RNASeq_Lasorsa/006_RNASeq_patients/006_res_shrunk.rds")
best_diff <- res_df %>% 
  filter(abs(log2FoldChange) > 5 | padj < 1e-150) %>% 
  mutate(stat = sign(log2FoldChange)*-log10(padj)) %>% 
  arrange(-stat) %>% 
  slice(-(26:(n()-25))) 

auto_tpm_orig <- read_table("data/RNASeq_Lasorsa/BGI_F21FTSEUHT2056_HOMypeiR/Expression/gene_expression.xls")
auto_tpm <- auto_tpm_orig %>% 
  select(gene_symbol, starts_with("tpm"), -contains("KO")) 

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

png("plots/006c/001_Heatmap_best_genes.png", h = 3000, w = 4200, res = 400)
Heatmap(best_diff_tpm,
        cluster_columns = F,
        col = colors_hm,
        column_labels = str_replace_all(colnames(best_diff_tpm), "_", " "),
        row_names_gp = grid::gpar(fontsize = 8),
        column_names_rot = 30,
        column_split = rep(1:2, each = 9),
        column_title = "Top and bottom best differentially expressed genes",
        heatmap_legend_param = list(
          title = "log10-scaled \nTPM")
        )
dev.off()

# Pathways
pathways <- read_delim("results/pathways_to_search_Lasorsa.csv", col_names = c("name", "x2", "gene_symbol")) %>% 
  select(1, 3)

# Remove TESSA samples and compute means
path_mean <- auto_tpm |> 
  filter(gene_symbol %in% pathways$gene_symbol) |> 
  rename_with(~ gsub("tpm_", "", .x)) %>%
  mutate(across(where(is.numeric), ~ log10(.x + 1))) |> 
  rowwise() |> 
  mutate(C1 = mean(c_across(starts_with("ctrATCC")), na.rm = T),
         C2 = mean(c_across(starts_with("ctrMS")), na.rm = T),
         P1 = mean(c_across(starts_with("agc_1")), na.rm = T),
         P2A = mean(c_across(starts_with("agc_2A")), na.rm = T),
         P2B = mean(c_across(starts_with("agc_2B")), na.rm = T)) %>%
  ungroup() |> 
  # filter(!if_all(where(is.numeric), ~ .x == 0)) %>%
  select(gene_symbol, C1, C2, P1, P2A, P2B)

pathway_un <- unique(pathways$name)

for(i in pathway_un){
  
  path_genes <- pathways %>% 
    filter(name == i) %>% 
    pull(gene_symbol)
  
  path_mtx <- 
    # auto_tpm %>%
    path_mean |> 
    filter(gene_symbol %in% path_genes) %>% 
    #arrange(as.numeric(str_extract(gene_symbol, "[0-9]{1,2}$"))) %>% 
    # mutate(across(where(is.numeric), ~ log10(.x + 1))) %>% 
    # filter(!if_all(where(is.numeric), ~ .x == 0)) %>% 
    # rename_with(~ gsub("tpm_", "", .x)) %>% 
    # mutate(C1 = across(starts_with("ctrATCC"), ~mean(.x)),
    #        C2 = across(starts_with("ctrATCC"), ~mean(.x)),)
    # relocate(starts_with("ctr")) %>% 
    column_to_rownames("gene_symbol") %>% 
    as.matrix() %>% 
    t() %>% 
    scale() %>% 
    t() |> 
    na.omit()
  
  if(nrow(path_mtx) > 30){
    row_fontsize <- 7
    plot_h <- 3000
    plot_w <- 3000
  }else{row_fontsize <- 10
        plot_h <- 1500
        plot_w <- 2500}
  
  colors_hm <- circlize::colorRamp2(
    breaks = seq(-1, 1, length = 3), 
    colors = c("red", "white", "blue"), 
    space = "RGB")
  
  pathway_file_name <- paste0(which(pathway_un == i), "_", str_replace_all(i, "[:punct:]|\\s", "_"))
  
  png(paste0("plots/006c/Heatmaps_pathway/", pathway_file_name, ".png"), h = plot_h, w = plot_w, res = 450)
  hm <- Heatmap(
    path_mtx,
    cluster_columns = F,
    column_names_centered = T,
    col = colors_hm,
    # column_labels = str_replace_all(colnames(best_diff_tpm), "_", " "),
    row_names_gp = grid::gpar(fontsize = row_fontsize),
    column_names_gp = grid::gpar(fontsize = 10),
    column_names_rot = 0,
    # column_split = rep(1:2, each = 9),
    column_split = c(1, 1, 2, 2, 2),
    column_title = paste0(i, " genes expression"),
    column_title_gp = grid::gpar(fontsize = 10),
    heatmap_legend_param = list(title = "log10-scaled \nTPM")
  )
  print(hm)
  dev.off()
}

#GSEA
library(fgsea)

mlist <- split(pathways$gene_symbol, pathways$name)
stat <- res_df %>% 
  mutate(stat = sign(log2FoldChange)*-log10(padj)) %>% 
  drop_na(stat) %>% 
  pull(stat, gene_symbol)

set.seed(41)
fgsea_res <- fgseaMultilevel(pathways = mlist,
                stats = stat,
                nproc = 6)
# The only significant pathway is the ATP synthase pathway
NES <- fgsea_res[pathway == "respiratory chain / ATP synthase"]$NES
padj <- fgsea_res[pathway == "respiratory chain / ATP synthase"]$padj

source("../../functions/plotEnrichment2.R")
png("plots/006c/GSEA_plot_ATP_synthase.png", h = 1500, w = 3000, res = 400)
plotEnrichment2(mlist[["respiratory chain / ATP synthase"]], stat)+
  labs(title = "respiratory chain/ATP synthase", 
       subtitle = paste0("NES = ", round(NES, 3), "  p.adj = ", scales::scientific(padj))) +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
        plot.subtitle = element_text(size = 8, hjust = 0.5)) +
  xlim(0,20000)
dev.off()

# Dotplots genes of interest ----
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

targets <- pathways %>% 
  pull(gene_symbol)

auto_new_l2fc <- auto_res_orig %>% 
  filter(gene_symbol %in% targets) %>%
  select(gene_symbol, contains("log2fc")) %>% 
  select(!contains("KO52")) %>% 
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
         padj = ifelse(padj == 0, .Machine$double.xmin, padj)) %>%
  ungroup()

all(rownames(auto_new_l2fc) == rownames(auto_new_padj)) # TRUE
all.equal(auto_new_l2fc$gene_symbol, auto_new_padj$gene_symbol) # TRUE

auto_new <- cbind(auto_new_l2fc, auto_new_padj %>% select(padj)) %>% 
  mutate(contrast = .$contrast %>% fct_relabel(~gsub("(MS|ATCC|TESSA|ctr)(-vs-)(agc|agc_1|agc_2A|agc_2B)", "\\3\\2\\1", .x)),
         contrast = as_factor(contrast) %>% 
           fct_relevel(., "agc-vs-ATCC", after = Inf) %>% 
           fct_relevel(., "agc-vs-MS", after = Inf) %>%
           fct_relevel(., "agc-vs-TESSA", after = Inf) %>% 
           fct_relevel(., "agc-vs-ctr", after = Inf) %>% 
           fct_relevel(., rev),
         contrast = fct_relabel(contrast, ~gsub("\\-", " ", .x))
  )

for(i in pathway_un){
  
  auto_target <- pathways %>% 
    filter(name == i) %>% 
    pull(gene_symbol)
  
  auto_new %>% 
    filter(gene_symbol %in% auto_target, !is.na(padj)) %>% 
    arrange(gene_symbol) %>% 
    mutate(logpadj = -log10(padj),
           logpadj = ifelse(logpadj >= 20, 20, logpadj)) %>% # Doing this because otherwise can't change the range of circles sizes using scale_size_continuous and large pvalues make all other supersmall
    ggplot(aes(gene_symbol, contrast, size = logpadj)) +
    geom_point(aes(fill = log2FC), shape = 21, color = "black") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = paste(i, " - Non aggregated and aggregated contrasts log2FoldChange"),
         subtitle = "Only points with padj < 0.05 are represented",
         x = "Gene",
         y = "Contrast") +
    scale_color_distiller(palette = "Spectral", limits = c(-2,2), oob = scales::squish, aesthetics = "fill") +
    scale_size_continuous(limits = c(-log10(0.05), 20), 
                          range = c(1, 5.5),
                          breaks = seq(0, 20, 5),
                          labels = c(as.character(seq(0, 15, 5)), "20+")
                          )
  
  ggsave(paste0("plots/006c/dotplots_expression/006c_log2FC_dotplot_", str_replace_all(i, "[:punct:]|\\s", "_"), ".png"), 
         h = 2500, w = 5200, dpi = 500, units = "px")
}
