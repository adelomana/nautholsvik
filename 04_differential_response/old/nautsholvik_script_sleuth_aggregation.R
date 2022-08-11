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
# 
# # install sleuth
# BiocManager::install("rhdf5")
# BiocManager::install("devtools")
# devtools::install_github("pachterlab/sleuth")

library(biomaRt)
library(sleuth)
library(crayon) # so the messages are blue
library(tictoc)

#
# 0. user-defined variables
#
setwd("~/scratch/")

kallisto_dir = "/home/adrian/projects/nautholsvik/results/kallisto/kallisto.100"
metadata_file = "/home/adrian/projects/nautholsvik/metadata/nautsholvik_metadata - Sheet1.tsv"
results_dir = '/home/adrian/projects/nautholsvik/results/DEGs_sleuth'

#
# 1. generate gene to transcript mapping
#

# annotation defined from sleuth walk through, https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'https://uswest.ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
View(t2g)

#
# 2. read metadata
#
metadata = read.table(metadata_file, sep='\t', header=TRUE)
View(metadata)

#
# 3. work on different contrasts
#
cat(blue('quantify the effect of siMITF, IFN and their interaction'), fill=TRUE)

# 3.1. prepare metadata including paths
rules = metadata$genotype != 'wt'
s2c = metadata[rules, ]

paths = file.path(kallisto_dir, s2c$sample, "abundance.h5")
s2c$path = paths
View(s2c)

# 3.2. prepare object for sleuth and plotting
cat(blue('preparing sleuth object...'), fill=TRUE)

tic()
so = sleuth_prep(s2c, target_mapping=t2g, aggregation_column='ens_gene') 
toc()

tic()
sop = sleuth_prep(s2c, 
                  target_mapping=t2g, 
                  aggregation_column='ens_gene', 
                  read_bootstrap_tpm=TRUE, 
                  extra_bootstrap_summary=TRUE, 
                  read_bootstrap_tpm=TRUE, 
                  gene_mode=TRUE)
toc()

# 3.3. fit all models 
cat(blue('building models...'), fill=TRUE)
tic()
so = sleuth_fit(so, ~genotype + treatment + genotype:treatment, 'interaction')
so = sleuth_fit(so, ~genotype + treatment, 'full')
so = sleuth_fit(so, ~genotype, 'treatment')
so = sleuth_fit(so, ~treatment, 'genotype')
#so = sleuth_fit(so, ~1, 'reduced')
toc()

# 3.4.1. statistical test about the effect of IFN
cat(blue('LRT testing...'), fill=TRUE)
lrt = sleuth_lrt(so, 'treatment', 'full')
lrt_table = sleuth_results(lrt, 'treatment:full', 'lrt', show_all=FALSE, pval_aggregate=TRUE)
lrt_table = dplyr::filter(lrt_table, qval <= 0.05)
print(dim(lrt_table))
head(lrt_table, 10)

plot_bootstrap(sop, lrt_table$target_id[1], units="tpm", color_by = "treatment")
plot_bootstrap(sop, lrt_table$target_id[2], units="tpm", color_by = "treatment")
plot_bootstrap(sop, lrt_table$target_id[3], units="tpm", color_by = "treatment")
plot_bootstrap(sop, lrt_table$target_id[dim(lrt_table)[1]/2], units="tpm", color_by = "treatment")
plot_bootstrap(sop, lrt_table$target_id[dim(lrt_table)[1]], units="tpm", color_by = "treatment")

cat(blue('storing...'), fill=TRUE)
write.csv(lrt_table, file.path(results_dir, paste('IFN_response', 'csv', sep='.')))

# 3.4.2. determine genes responding to siMITF
cat(blue('LRT testing...'), fill=TRUE)
lrt = sleuth_lrt(so, 'genotype', 'full')
lrt_table = sleuth_results(lrt, 'genotype:full', 'lrt', show_all=FALSE, pval_aggregate=TRUE)
lrt_table = dplyr::filter(lrt_table, qval <= 0.05)
print(dim(lrt_table))
head(lrt_table, 10)

plot_bootstrap(sop, lrt_table$target_id[1], units="tpm", color_by = "genotype")
plot_bootstrap(sop, lrt_table$target_id[2], units="tpm", color_by = "genotype")
plot_bootstrap(sop, lrt_table$target_id[3], units="tpm", color_by = "genotype")
plot_bootstrap(sop, lrt_table$target_id[dim(lrt_table)[1]/2], units="tpm", color_by = "genotype")
plot_bootstrap(sop, lrt_table$target_id[dim(lrt_table)[1]], units="tpm", color_by = "genotype")
plot_bootstrap(sop, 'ENSG00000187098', units="tpm", color_by = "genotype") # MITF

cat(blue('storing...'), fill=TRUE)
write.csv(lrt_table, file.path(results_dir, paste('siMITF_response', 'csv', sep='.')))

# 3.4.3. determine genes responding non-linearly
cat(blue('LRT testing...'), fill=TRUE)
lrt = sleuth_lrt(so, 'full', 'interaction')
lrt_table = sleuth_results(lrt, 'full:interaction', 'lrt', show_all=FALSE, pval_aggregate=TRUE)
lrt_table = dplyr::filter(lrt_table, qval <= 0.05)
print(dim(lrt_table))
head(lrt_table, 10)

plot_bootstrap(sop, lrt_table$target_id[1], units="tpm", color_by = "genotype")
plot_bootstrap(sop, lrt_table$target_id[2], units="tpm", color_by = "genotype")
plot_bootstrap(sop, lrt_table$target_id[3], units="tpm", color_by = "genotype")
plot_bootstrap(sop, lrt_table$target_id[dim(lrt_table)[1]/2], units="tpm", color_by = "genotype")
plot_bootstrap(sop, lrt_table$target_id[dim(lrt_table)[1]], units="tpm", color_by = "genotype")
plot_bootstrap(sop, 'ENSG00000120217', units="tpm", color_by = "genotype") # PD-L1

cat(blue('storing...'), fill=TRUE)
write.csv(lrt_table, file.path(results_dir, paste('interaction_response', 'csv', sep='.')))



