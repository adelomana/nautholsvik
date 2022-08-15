library(biomaRt)

#
# 0. user-defined variables
#
setwd("~/scratch/")

annotation_file = "/home/adrian/projects/nautholsvik/results/annotation/annotation.csv"

#
# 1. generate gene to transcript mapping
#

# annotation defined from sleuth walk through, https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
mart = biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host='https://uswest.ensembl.org')
t2g = biomaRt::getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart=mart)
t2g = dplyr::rename(t2g, target_id=ensembl_transcript_id, ens_gene=ensembl_gene_id, ext_gene=external_gene_name)

write.csv(t2g, annotation_file)