library(tidyverse)
library(Rsubread)
library(Rsamtools)

# Align RNASeq FASTQ using Rsubread
genome <- "data/ATACSeq/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
rsub_ind <- gsub("\\.fa.gz$", "", genome)

buildindex(rsub_ind, genome)

fq <- list.files("../agc1/data/delivery_20200428/raw_reads/", 
                 pattern = "fastq.gz", full.names = TRUE)

for(rname in fq){
  
  outBAM <- paste0("data/007/bams_rnaseq/", 
                   basename(rname) %>% str_extract("Olineu[0-9]{1}"),
                   ".bam")
  align(rsub_ind,
        readfile1 = rname,
        readfile2 = gsub("_R1_", "_R2_", rname),
        output_file = outBAM,
        nthreads = 5,
        type = "rna",
        unique = TRUE)
  
  sortedBAM <- str_replace(outBAM, "\\.bam", ".sorted.bam")
  sortBam(outBAM, gsub("\\.bam", "", sortedBAM))
  
  indexBam(sortedBAM)
}

# Get coverage + GC content and plots ----
# Run 007_bamCoverge.sh
files <- list.files("data/ATACSeq/new_bam", pattern = "*.bedgraph", full.names = TRUE)

for(file in files){

  bed <- read_delim(file, col_names = c("Chromosome", "Start", "End", "Count"))
  sample <- str_extract(file, "_[:alnum:]{2}_[:alnum:]{3}.") %>%
    str_replace_all("[:punct:]", " ")

  bedfil <- bed %>%
    mutate(Region = rowMeans(select(bed, Start, End), na.rm = T), .after = Chromosome) %>%
    filter(Chromosome %in% c(seq(1:19), "X"), Count != 0) %>%
    rowwise() %>%
    mutate(Chromosome = paste0("chr", Chromosome),
           Chromosome = factor(Chromosome, levels = paste0("chr", c(seq(1:19), "X"))),
           count_norm = (Count - mean(Count))/sd(Count))

  bedgc <- bedfil %>%
    bind_cols(gc) %>%
    rename(GC_content = last_col()) %>%
    mutate(count_zscore = (Count - count_me)/count_sd,
           gc_zscore = (GC_content - gc_me)/gc_sd) %>% # Compute Z score of Counts to remove outliers
    filter(abs(count_zscore) < 3, abs(gc_zscore) < 3) %>%
    group_by(Chromosome) %>%
    mutate(p_val = cor.test(log10(Count), GC_content)$p.value,
           p_val = replace(p_val, p_val == 0e+00, 0 + .Machine$double.eps),
           p_adj = p.adjust(p_val, method = "bonferroni") %>% scales::scientific(),
           adj_rsq = summary(lm(Count ~ GC_content))$adj.r.squared %>% signif(3),
           chr_label = paste0(Chromosome, " - ",
                          "p_adj = ", p_adj, " - ",
                          "adj_rsq = ", adj_rsq)) %>%
    ungroup()

  ggplot(bedgc, aes(x = Region)) +
    geom_point(aes(y = Count/mean(Count)), size = 1, color = "black", alpha = 0.2) +
    geom_point(aes(y = GC_content/mean(GC_content)), size = 1, color = "#ED254E", alpha = 0.1)  +
    scale_y_continuous(name = "Normalized Counts",
                       sec.axis = sec_axis(trans = ~., name = "Normalized GC content"),
                       trans = "log10") +
    facet_wrap(~Chromosome, scales = "free") +
    theme_light() +
    labs(title = paste0("Counts and GC content distribution in sample", sample)) +
    theme(plot.title = element_text(size = 25),
          axis.title.x = element_text(size = 18),
          axis.title.y.left = element_text(color = "black", size = 18),
          axis.text.y.left = element_text(color = "black"),
          axis.title.y.right = element_text(color = "#ED254E", size = 18),
          axis.text.y.right = element_text(color = "#ED254E"),
          strip.text.x = element_text(color = "black", size = 12),
          strip.background = element_rect(fill = "#ED254E")
    )

  ggsave(paste0("plots/005_ATACSeq/005.1_3_Counts_and_GC_content_over_region",
                str_replace_all(sample, "\\s", "_"),
                ".png"),
         h = 35,
         w = 70,
         units = "cm",
         dpi = 600)

  ggplot(bedgc, aes(GC_content, Count)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#4775FF") +
    scale_y_log10() +
    facet_wrap(~chr_label, scales = "free") +
    theme_light() +
    labs(title = paste0("Counts distribution over GC content in sample", sample)) +
    theme(plot.title = element_text(size = 25),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          strip.background = element_rect(fill = "#4775FF"),
          strip.text = element_text(size = 12))


    ggsave(paste0("plots/005_ATACSeq/005.1_4_Counts_over_GC_content",
                  str_replace_all(sample, "\\s", "_"),
                  ".png"),
           h = 35,
           w = 70,
           units = "cm",
           dpi = 600)
}