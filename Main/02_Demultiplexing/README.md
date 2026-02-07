# Demultiplexing and doublet detection

To assign the nuclei to each donor and identify doublets, we used Vireo implemented with Demuxafy (version: 2.1.0)

We first used `generate_sample_list.R` and `subset_pool_vcf.qsub.sh` to subset the genotype VCF files for each pool

Then we used `vireo_TOB_ATAC_pool_1.qsub.sh` and `vireo_TOB_ATAC_pool_2.qsub.sh` to run Viero for each of the repeated pools. The input files are directly from the Cell Ranger ATAC output, so they should be reproducible.  

The output files will be used in the next step to remove doublets and unassigned cells/nuclei

