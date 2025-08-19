#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# args <- c("CD4_TCM", "21")
cell_type <- args[1]
chrNumber <- args[2]

# Upload libraries
library(data.table)
library(dplyr)
library(rtracklayer)
library(stringr)
setDTthreads(4)

# Directories
work_dir <- "/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/SMR/"

# cell_type_list <- sort(list.dirs(paste0(work_dir, "output/main_analyses/preprocessing_caQTL_rawID/"), full.names = FALSE, recursive = FALSE))

# for (cell_type in cell_type_list) {
#   
besd_dir <- paste0(work_dir, "output/main_analyses_NewGenotypes/preprocessing_caQTL_rawID/", cell_type, "/BESD/")

# Read in EPI
little_epi_filepath <- paste0(besd_dir, "Chr", chrNumber, ".epi")
# What this needs: chromosome, probeID, genetic distance, physical position, gene name, gene orientation
epi_df <- fread(little_epi_filepath)

# Output example from Anne
# seqid              V2    V3     start       gene_name strand
# <char>          <char> <int>     <int>          <char> <char>
#   1:      1 ENSG00000000457     0 169849631           SCYL3      -
#   2:      1 ENSG00000000460     0 169662007        C1orf112      +
#   3:      1 ENSG00000000938     0  27612064             FGR      -

# Start to reformulate the .epi file
setnames(epi_df, "V1", "seqid")
# Convert from <lgcl> to <char>
epi_df[, seqid := as.character(seqid)]
# Modify the specific value
epi_df[, seqid := as.character(chrNumber)]

# setnames(epi_df, "V4", "start")
# epi_df[, start := as.integer(start)]
# epi_df[, start := sub(".*:(\\d+)-\\d+", "\\1", V2)]
split_colon <- str_split(epi_df$V2, ":", simplify = TRUE)
positions <- str_split(split_colon[, 2], "-", simplify = TRUE)
position1 <- as.numeric(positions[, 1])
position2 <- as.numeric(positions[, 2])
epi_df$V4 <- as.integer(round((position1 + position2) / 2))

epi_df$V5 = epi_df$V2
epi_df[, V6 := as.character(V6)]
epi_df[, V6 := "NA"]

fwrite(epi_df, little_epi_filepath, sep="\t", quote=F, col.names=F)

# Our output would look like
# seqid                       V2    V3     start
# <char>                   <char> <int>    <char>
# 1:     10       chr1:816955-817603     0    816955
# 2:     10       chr1:821186-821701     0    821186
# 3:     10       chr1:841229-841821     0    841229

# }


#for (cell_type in cell_type_list) {
  
# besd_dir <- paste0(work_dir, "output/main_analyses/preprocessing_caQTL_rawID/", cell_type, "/BESD/")
# for (chrNumber in seq(1,22)){
#   # Read in EPI
#   little_epi_filepath <- paste0(besd_dir, "Chr", chrNumber, ".epi")
#   if(!file.exists(little_epi_filepath)) {
#     next  # Skip to the next chromatin
#   }
#   # What this needs: chromosome, probeID, genetic distance, physical position, gene name, gene orientation
#   epi_df <- fread(little_epi_filepath)
#   
#   split_colon <- str_split(epi_df$V2, ":", simplify = TRUE)
#   positions <- str_split(split_colon[, 2], "-", simplify = TRUE)
#   position1 <- as.numeric(positions[, 1])
#   position2 <- as.numeric(positions[, 2])
#   epi_df$V4 <- as.integer(round((position1 + position2) / 2))
#   
#   fwrite(epi_df, little_epi_filepath, sep="\t", quote=F, col.names=F)
# }
# #}