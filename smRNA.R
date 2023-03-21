library(ggpubr)
library(tidyverse)
library(ggplot2)
library(lemon)
library(scales)
library(RColorBrewer)
library(ggrepel)
library(ggpmisc)
library(ggseqlogo)
library(rstatix)
library(ggdist)

source("/fs/ess/PCON0160/ben/bin/mighty.R")

##########################################
# plot 22G level eggd-1 vs. wild type
eggd1_dat = read.delim("smRNA_results/DGE/rfppgl_vs_eggd1_total_norm.tsv")

p = xy_dge(eggd1_dat %>% filter(feature == "siRNA"), "rfppgl", "eggd1", 2^-10, 2^15, fold_change = 2); p 
p
ggsave(p, filename = "figs/siRNA_eggd1_vs_wt.png", dpi = 300, height = 10, width = 10, units = "cm")

##########################################
# plot piRNA level eggd-1 vs. wild type
eggd1_dat = read.delim("smRNA_results/DGE/rfppgl_vs_eggd1_total_norm.tsv")

p = xy_dge(eggd1_dat %>% filter(feature == "piRNA"), "rfppgl", "eggd1", 2^-10, 2^15, fold_change = 1); p 

eggd1_dat_long = eggd1_dat %>% 
  select(eggd1, rfppgl, feature, locus_id, class) %>% 
  pivot_longer(!c(feature, locus_id, class), names_to = 'sample', values_to = 'count') %>% 
  mutate(class = ifelse(grepl("Sperm", class), "Spermatogenic", 
                        ifelse(grepl("Oogeni", class), "Oogenic", "Other")))

eggd1_dat_long$sample_f = factor(eggd1_dat_long$sample, levels = c("rfppgl", "eggd1"))

p = plot_boxplot(eggd1_dat_long %>% filter(feature == "piRNA"), 
                 counts_col = 'count', samples_col = 'sample_f', 
                 ylog2 = T, distribution = T, pvals = T) + 
  theme(aspect.ratio = 1.5)
ggsave(p, filename = "figs/piRNA_eggd1_vs_wt_boxplot.png", dpi = 300, height = 4.5, width = 4.5)

# plot piRNA level for different subsets of piRNAs expression by sex
p = plot_boxplot(eggd1_dat_long %>% filter(feature == "piRNA"), 
                 level_col = 'class',
                 counts_col = 'count', 
                 samples_col = 'sample_f', 
                 ylog2 = T, distribution = T, pvals = T) + 
  theme(aspect.ratio = 1.5)
p
ggsave(p, filename = "figs/piRNA_eggd1_vs_wt_boxplot_sex.png", dpi = 300, height = 4.5, width = 10)


# plot piRNA ox vs. unox
counts = read.delim("smRNA_results/master_tables/eggd_analysis.aligned.v0.m1000.count_total_norm.average.tsv")

counts_long = counts %>% 
  select(-biotype, -gene_name, -class, -seq_id) %>% 
  pivot_longer(!c(feature, locus_id), names_to = 'sample', values_to = 'count')
levels(factor(counts_long$sample))

conditions = c("N2Oxi", "rfppgl", "eggd1", "eggd1Oxi")

p1 = plot_boxplot(counts_long %>% filter(feature == "piRNA") %>% filter(sample %in% c("N2Oxi", "rfppgl")), 
                 counts_col = 'count', samples_col = 'sample', 
                 ylog2 = T, distribution = T, pvals = T) + 
  theme(aspect.ratio = 1.5)

p2 = plot_boxplot(counts_long %>% filter(feature == "piRNA") %>% filter(sample %in% c("eggd1", "eggd1Oxi")), 
                  counts_col = 'count', samples_col = 'sample', 
                  ylog2 = T, distribution = T, pvals = T) + 
  theme(aspect.ratio = 1.5)

p = ggarrange(p1, p2, ncol = 2)
ggsave(p, filename = "figs/piRNA_eggd1_vs_wt_boxplot_oxidation.png", dpi = 300, height = 8, width = 8)


##########################################
# plot alg/ergo target 26G level
counts = read.delim("smRNA_results/master_tables/eggd_analysis.aligned.v0.m1000.count_total_norm.average.tsv")

counts_long = counts %>% 
  select(-biotype, -gene_name, -seq_id) %>% 
  pivot_longer(!c(feature, locus_id, class), names_to = 'sample', values_to = 'count')

levels(factor(counts_long$sample))

conditions = c("rfppgl", "eggd1")

p1 = plot_boxplot(counts_long %>% filter(feature == "siRNA26G" & grepl("ALG", class)),
                  conditions = conditions,
                  counts_col = 'count', 
                  samples_col = 'sample', 
                  ylog2 = T, ymin = 2^-10, ymax = 2^15,
                  distribution = T, 
                  pvals = T) + 
  theme(aspect.ratio = 1.5) 
p1

p2 = plot_boxplot(counts_long %>% filter(feature == "siRNA26G" & grepl("ERGO", class)),
                  conditions = conditions,
                  counts_col = 'count', 
                  samples_col = 'sample', 
                  ylog2 = T, ymin = 2^-10, ymax = 2^15, dots = T,
                  distribution = T, 
                  pvals = T) + 
  theme(aspect.ratio = 1.5) 
p2

p = ggarrange(p1, p2, ncol = 2)
ggsave(p, filename = "figs/26G_level_eggd1_vs_wt_boxplot.png", dpi = 300, height = 8, width = 8)


##########################################
# length distribution wild type vs eggd 1
replicates = read.delim("smRNA_results/samples/replicates.csv", sep = ',')

lendist = read.delim("smRNA_results/eggd_lendis.tsv")

lendist$first_nt = as.character(lendist$first_nt)

lendist = lendist %>% 
  mutate(first_nt = ifelse(first_nt == TRUE, "T", first_nt))

lendist_grouped = lendist %>% 
  left_join(replicates, by = c("sample" = "simple_name")) %>% 
  group_by(seqlen, first_nt, condition) %>% 
  summarise(M = mean(sum)) %>% 
  group_by(condition) %>% 
  mutate(percent = 100 * (M / sum(M) ))

sub_lendis_grouped = lendist_grouped %>% filter(condition == "eggd1" || condition == "rfppgl")

sub_lendis_grouped$condition_f = factor(sub_lendis_grouped$condition, levels = c("rfppgl", "eggd1"), ordered = TRUE)

p = ggplot(data = sub_lendis_grouped, aes(x = seqlen, y = M, fill = condition_f)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  my_theme() + 
  theme(aspect.ratio = 0.5) + 
  scale_x_continuous(limits = c(15,30), breaks = c(seq(15,30,by=3))) + 
  scale_fill_manual(values = c("rfppgl" = "blue", "eggd1" = "red"))
p
ggsave(p, filename = "figs/piRNA_eggd1_vs_wt_length_distribution.png", dpi = 300, height = 4.5, width = 4.5)

##########################################
# plot the overlap of up/down siRNAs with gene classes

overlap_signif = function(N_genes_total, N_class_total, N_genes_subset, N_class_subset){
  
  p = fisher.test(matrix(c(N_class_subset,
                       N_genes_subset - N_class_subset, 
                       N_class_total - N_class_subset,
                       N_genes_total - N_genes_subset - N_class_total + N_class_subset), nrow = 2))$p.value
  
  i_total = N_class_total / N_genes_total
  i_sub = N_class_subset / N_genes_subset
  
  rep = log2( round( i_sub / i_total, 2) )
  
  return(c(rep, p))
  
}

siRNA_classes = c("HRDE1", "CSR1", "PRG1", "WAGO1", "MUT16", "WAGO4")
pathways_all = read.delim("/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/pathways.tsv", sep = "\t", col.names = c("gene_name", "locus_id", "seq_id", "biotype", "class"))

eggd1_dat = read.delim("smRNA_results/DGE/rfppgl_vs_eggd1_total_norm.tsv")

down = eggd1_dat %>% filter(feature == "siRNA") %>% filter(log2FoldChange <= -1 & pvalue < 0.05)
up =  eggd1_dat %>% filter(feature == "siRNA") %>% filter(log2FoldChange >= 1 & pvalue < 0.05)


total_genes = c()
total_in_class = c()
total_down = c()
down_in_class = c()
down_enrichment = c()
down_pvalue = c()
total_up = c()
up_in_class = c()
class_name = c()
up_enrichment = c()
up_pvalue = c()
i = 1
for (c in siRNA_classes){
  
  # total genes in genome
  N_genes = nrow(pathways_all)
  
  # total genes in class 
  N_class = pathways_all %>% filter(grepl(c, class)) %>% nrow()
  
  # total genes down
  N_down_all = down %>% nrow()
  
  # total genes down in class
  N_down_class = down %>% filter(grepl(c, class)) %>% nrow()
  
  # total genes up
  N_up_all = up %>% nrow()
  
  # total genes up in class 
  N_up_class = up %>% filter(grepl(c, class)) %>% nrow()
  
  # up stats
  up_stats = overlap_signif(N_genes, N_class, N_up_all, N_up_class)
  up_enrich = up_stats[[1]]
  up_p = up_stats[[2]]
  
  # down stats
  down_stats = overlap_signif(N_genes, N_class, N_down_all, N_down_class)
  down_enrich = down_stats[[1]]
  down_p = down_stats[[2]]  
  
  class_name[i] = c
  total_genes[i] = N_genes
  total_in_class[i] = N_class
  
  down_in_class[i] = N_down_class
  down_enrichment[i] = down_enrich
  down_pvalue[i] = down_p
  
  up_in_class[i] = N_up_class
  up_enrichment[i] = up_enrich
  up_pvalue[i] = up_p
  
  i = i + 1
  
}

pathway_counts = data.frame(
  "class_name" = class_name,
  "total_genes" = total_genes,
  "total_in_class" = total_in_class,
  "total_down" = N_down_all,
  "down_in_class" = down_in_class,
  "down_enrichment" = down_enrichment,
  "down_pvalue" = down_pvalue,
  "total_up" = N_up_all,
  "up_in_class" = up_in_class,
  "up_enrichment" = up_enrichment,
  "up_pvalue" = up_pvalue
)

add_psignif = function(pval){
  
  if (pval >= 0.05){
    return("ns")
  } else if (pval < 0.0001) {
    return("****")
  } else if (pval < 0.001) {
    return("***")
  } else if (pval < 0.01){
    return("**")
  } else {
    return("*")
  }
  
}

pathway_counts

pathway_counts_stacked = pathway_counts %>% 
  select(class_name, down_enrichment, down_pvalue, up_enrichment, up_pvalue) %>% 
  pivot_longer(!class_name, names_to = "sample", values_to = "val") %>% 
  mutate(group = ifelse(grepl("down", sample), "down", "up")) %>% 
  mutate(info = ifelse(grepl("enrich", sample), "fold_enrichment", "pvalue")) %>% 
  select(-sample) %>% 
  pivot_wider(names_from = info, values_from = val) %>% 
  mutate(compartment = ifelse(grepl("HRDE", class_name), "Nuclear", 
                              ifelse(grepl("CSR1|PRG1|WAGO1", class_name), "P_granule",
                                     ifelse(grepl("WAGO4", class_name), "Z_granule", "M_focus")))) %>% 
  rowwise() %>% 
  mutate(psignif = add_psignif(pvalue))

fills = c(
  "P_granule" = "springgreen3",
  "M_focus" = "mediumorchid1",
  "Nuclear" = "orange",
  "Z_granule" = "cyan"
)

pathway_counts_stacked$compartment_f = factor(pathway_counts_stacked$compartment, levels = c("Nuclear", "P_granule", "Z_granule", "M_focus"))

p_up = ggplot(data = pathway_counts_stacked %>% filter(group == "up"), aes(x = class_name, y = fold_enrichment, fill = compartment_f)) +
  geom_bar(color = "black", stat = "identity", width = 0.6) +
  my_theme_free_aspect()+
  scale_y_continuous(limits = c(-3,3), breaks = seq(-3,3,by=1)) + 
  #scale_x_discrete(expand = c(.5,.5)) + 
  scale_fill_manual(values = fills) + 
  scale_color_manual(values = fills) + 
  facet_grid(.~compartment_f, space = "free", scales = "free") +
  geom_text(aes(label=psignif, y = ifelse(fold_enrichment>0, fold_enrichment+0.1, fold_enrichment-0.2)), vjust=0) + 
  geom_hline(yintercept = 0, color = "grey60", lty = "dashed") + 
  ggtitle("Up-regulated")

p_down = ggplot(data = pathway_counts_stacked %>% filter(group == "down"), aes(x = class_name, y = fold_enrichment, fill = compartment_f)) +
  geom_bar(color = "black", stat = "identity", width = 0.6) +
  my_theme_free_aspect()+
  scale_y_continuous(limits = c(-3,3), breaks = seq(-3,3,by=1)) + 
  scale_fill_manual(values = fills) + 
  scale_color_manual(values = fills) + 
  facet_grid(.~compartment_f, space = "free", scales = "free") + 
  geom_text(aes(label=psignif, y = ifelse(fold_enrichment>0, fold_enrichment+0.1, fold_enrichment-0.2)), vjust=0) + 
  geom_hline(yintercept = 0, color = "grey60", lty = "dashed") + 
  ggtitle("Down-regulated")


p = ggarrange(p_up, p_down, align = 'v', ncol = 1)
p
ggsave(p, filename = "22G_de_overlap_pathways_4fold.pdf", dpi = 300, height = 6, width = 15)






