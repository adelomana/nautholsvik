#
# -1. packages installation
#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# setRepositories(ind=c(1:6))
# 
# BiocManager::install("biomaRt")

library(biomaRt)
library(tximport)
library(DESeq2)
library(BiocParallel)
library(crayon) # so the messages are blue

#
# 0. user-defined variables
#
register(MulticoreParam(20))
setwd("~/scratch/")

kallisto_dir = "/home/adrian/projects/nautholsvik/results/kallisto/kallisto.100"
metadata_file = "/home/adrian/projects/nautholsvik/metadata/nautsholvik_metadata - Sheet1.tsv"
results_dir = '/home/adrian/projects/nautholsvik/results/DEGs_DESeq2'

threshold = 10

#
# 1. generate gene to transcript mapping
#
# annotation defined from sleuth walk through, https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
mart = biomaRt::useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", version=96)
t2g = biomaRt::getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart=mart)
t2g = dplyr::rename(t2g, target_id=ensembl_transcript_id, ens_gene=ensembl_gene_id, ext_gene=external_gene_name)

#
# 2. read metadata
#
metadata = read.table(metadata_file, sep='\t', header=TRUE)
View(metadata)

#
# 3. prepare metadata including paths and read data
#
rules = metadata$genotype != 'wt'
working_metadata = metadata[rules, ]
View(working_metadata)

files = file.path(kallisto_dir, working_metadata$sample, "abundance.h5")
cat(blue('files'), fill=TRUE)
print(files)

txi = tximport(files, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

#
# 4. model comparisons
#

# #
# # 4.1. call significance on influences due to genotype only
# #
# dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~genotype) 
# dds$genotype = relevel(dds$genotype, ref="control")
# 
# keep = rowMaxs(counts(dds)) >= threshold
# dds = dds[keep, ]
# 
# dds = DESeq(dds, test="LRT", reduced=~1)
# 
# res = results(dds, parallel=TRUE) 
# filt1 = res[which(res$pvalue < 0.05), ]
# filt2 = filt1[which(filt1$padj < 0.1), ]
# filt3 = filt2[which(abs(filt2$log2FoldChange) > 1), ]
# cat(blue(paste('DEGs found for genotype only:', dim(filt3)[1], sep=' ')), fill=TRUE)
# write.table(filt3, file=paste(results_dir, '/genotype_only.tsv', sep=''), quote=FALSE, sep='\t')
# 
# #
# # 4.2. call significance on influences due to treatment only
# #
# dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~treatment) 
# dds$treatment = relevel(dds$treatment, ref='control')
# 
# keep = rowMaxs(counts(dds)) >= threshold
# dds = dds[keep, ]
# 
# dds = DESeq(dds, test="LRT", reduced=~1)
# 
# res = results(dds, parallel=TRUE) 
# filt1 = res[which(res$pvalue < 0.05), ]
# filt2 = filt1[which(filt1$padj < 0.1), ]
# filt3 = filt2[which(abs(filt2$log2FoldChange) > 1), ]
# cat(blue(paste('DEGs found for treatment only', dim(filt3)[1], sep=' ')), fill=TRUE)
# write.table(filt3, file=paste(results_dir, '/treatment_only.tsv', sep=''), quote=FALSE, sep='\t')
# 
# #
# # 4.3. call significance on two factors better than one
# # 
# dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~genotype + treatment) 
# dds$genotype = relevel(dds$genotype, ref="control")
# dds$treatment = relevel(dds$treatment, ref='control')
# 
# keep = rowMaxs(counts(dds)) >= threshold
# dds = dds[keep, ]
# 
# dds = DESeq(dds, test="LRT", reduced=~genotype)
# 
# res = results(dds, parallel=TRUE) 
# filt1 = res[which(res$pvalue < 0.05), ]
# filt2 = filt1[which(filt1$padj < 0.1), ]
# filt3 = filt2[which(abs(filt2$log2FoldChange) > 1), ]
# cat(blue(paste('DEGs found for two factors better than genotype', dim(filt3)[1], sep=' ')), fill=TRUE)
# write.table(filt3, file=paste(results_dir, '/two_better_than_genotype.tsv', sep=''), quote=FALSE, sep='\t')
# 
# ###
# 
# dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~treatment + genotype) 
# dds$genotype = relevel(dds$genotype, ref="control")
# dds$treatment = relevel(dds$treatment, ref='control')
# 
# keep = rowMaxs(counts(dds)) >= threshold
# dds = dds[keep, ]
# 
# dds = DESeq(dds, test="LRT", reduced=~treatment)
# 
# res = results(dds, parallel=TRUE) 
# filt1 = res[which(res$pvalue < 0.05), ]
# filt2 = filt1[which(filt1$padj < 0.1), ]
# filt3 = filt2[which(abs(filt2$log2FoldChange) > 1), ]
# cat(blue(paste('DEGs found for two factors better than treatment', dim(filt3)[1], sep=' ')), fill=TRUE)
# write.table(filt3, file=paste(results_dir, '/two_better_than_treatment.tsv', sep=''), quote=FALSE, sep='\t')

#
# 4.4. call significance on the interaction: no fold-change filter
#
dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~genotype + treatment + genotype:treatment) 
dds$genotype = relevel(dds$genotype, ref="control")
dds$treatment = relevel(dds$treatment, ref='control')

keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]

dds = DESeq(dds, test="LRT", reduced=~genotype + treatment)

res = results(dds, parallel=TRUE) 
filt1 = res[which(res$pvalue < 0.05), ]
filt2 = filt1[which(filt1$padj < 0.1), ]
cat(blue(paste('DEGs found', dim(filt2)[1], sep=' ')), fill=TRUE)
write.table(filt2, file=paste(results_dir, '/interaction.tsv', sep=''), quote=FALSE, sep='\t')




