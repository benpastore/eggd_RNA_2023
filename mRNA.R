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

# counts
counts = read.table("mRNA_analysis/counts/analysis.tsv", header = TRUE, fill = TRUE, check.names=FALSE)
counts = counts %>% filter(!grepl("snRNA|snoRNA|scRNA|tRNA|asRNA|miRNA|piRNA", biotype)) %>% filter(biotype != "ncRNA")
colnames(counts) = gsub(".tsv","",colnames(counts))
head(counts)

# pathways 
pathways = read.delim("/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/pathways.tsv", sep = "\t",
                      col.names = c("gene_name", "seq_id", "locus_id", "biotype", "class"))
pathways = pathways %>% select(-seq_id, -locus_id)

# merge counts to pathways
counts = counts %>% 
  select(-class, -biotype) %>% 
  left_join(pathways, by = c("gene_id" = "gene_name"))

# coldata
coldata = read.table("mRNA_analysis/samples/replicates.csv", header = TRUE, sep = ",", col.names = c("sample", "condition"))
row.names(coldata) = coldata$sample
coldata = subset(coldata, select = c("condition"))
coldata

# select cols that are non-numeric/numeric
tmp = data.frame(counts, row.names = 1, check.names=FALSE)
names_df = select_if(tmp, Negate(is.numeric))
names_df['gene_id'] = rownames(names_df)
cts = select_if(tmp, is.numeric)
cts = cts[ , sort(colnames(cts))]
head(cts)

# remove rows that have no counts
x = rowSums(cts) >= 1
cts_filt = cts[x,]
cts_filt[] <- sapply(cts_filt, as.integer)
list(colnames(cts_filt))[[1]]

# set condition equal to coldata condition
condition = coldata$condition

# make sure columns in counts are in same order as they appear in coldata
rownames(coldata)
cts_filt = cts_filt[,rownames(coldata)]
all( rownames(coldata) == colnames(cts_filt) )

# DESeq
deobj <- DESeqDataSetFromMatrix(countData = cts_filt, colData = coldata, design = ~condition)
dds <- DESeq(deobj)
tmp = counts(dds, normalized = TRUE)
dds_counts = as.data.frame(tmp) %>% rownames_to_column(var = "gene_id")

counts_table = dds_counts %>% 
  left_join(names_df, by = "gene_id") %>% 
  select(gene_id, seq_id, locus_id, biotype, class, everything())

write.table( counts_table, "normalized_counts.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# PCA plot to assess variance within sample groups and between sample groups 
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

PCA = ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(aspect.ratio = 1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
PCA

#results and plot for each comparison
comparisons = data.frame(x = c("control"), y = c("eggd1"))
head(comparisons)

for(i in 1:nrow(comparisons)){
  
  sample_x = as.character(comparisons[i,1])
  sample_y = as.character(comparisons[i, 2])
  
  curr = paste0(sample_x,"-vs-",sample_y)
  print(curr)
  
  xdf = subset(coldata, condition == paste0(sample_x))
  xnames = rownames(xdf)
  
  ydf = subset(coldata, condition == paste0(sample_y))
  ynames = rownames(ydf)
  
  res = results(dds, contrast = c("condition",paste0(sample_y),paste0(sample_x)))
  
  c = as.data.frame(counts(dds,normalized = TRUE))
  keep = c(paste0(xnames), paste0(ynames))
  
  c_sub = subset(c, select = keep)
  
  resdata = merge(as.data.frame(res), c_sub, by = 'row.names', sort = FALSE)
  names(resdata)[1] <- "gene_id"
  
  resdata[paste0(sample_x)] = rowMeans( resdata[ , xnames] )
  resdata[paste0(sample_y)] = rowMeans( resdata[ , ynames] )
  
  resdata = merge(resdata, names_df, by = "gene_id", all.x = TRUE)
  
  deseq_cols = c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  cols = c(names_df, deseq_cols)
  
  #ord = c(cols, deseq_cols, list(xnames)[[1]], list(ynames)[[1]], sample_x, sample_y)
  write.table(resdata, paste0(curr,"_deseq_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  p = xy_dge(resdata, paste0(sample_x), paste0(sample_y), 2^-5, 2^20, deseq = TRUE, fold_change = 2)
  ggsave(p, filename = paste0("figs/",curr,"_deseq_results.png"), dpi = 300, height = 10, width = 10, units = "cm")
  
}
p



