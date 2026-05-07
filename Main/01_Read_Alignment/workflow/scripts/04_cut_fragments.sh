#!/bin/bash
#PBS -P <HPC_PROJECT>
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l storage=gdata/<HPC_PROJECT>+gdata/<HPC_PROJECT_2>
#PBS -l mem=128GB
#PBS -l ncpus=32
#PBS -l jobfs=100GB
#PBS -l wd
#PBS -e /g/data/<HPC_PROJECT>/repos/tenk10k_cellranger/logs/pbs_logs/cutfrags.stderr
#PBS -o /g/data/<HPC_PROJECT>/repos/tenk10k_cellranger/logs/pbs_logs/cufrags.stdout
#PBS -M <YOUR_EMAIL>
#PBS -m ae


module load parallel

. /home/<HPC_GROUP>/<HPC_USERNAME>/micromamba/etc/profile.d/micromamba.sh
micromamba activate genomic-utils


# This script removes the last column from fragment files, creating a copy in the same directory.
# A newer version of cellranger updated the fragment file format with an extra strand column;
# this script removes it to make fragment files compatible with downstream Signac / SnapATAC analysis.
# UPDATE: newer. versions of Signac / SnapATAC can handle the new format so this is not necessary anymore 
parallel -j32 'if [ ! -f {//}/fragments_nostrand.tsv.gz ]; then zcat {} | cut -f1-5 | bgzip > {//}/fragments_nostrand.tsv.gz; fi' ::: /g/data/<HPC_PROJECT>/data/atac/atac_count_outs/force_cells_test/*/*/outs/fragments.tsv.gz
parallel -j32 'if [ ! -f {}.tbi ]; then tabix -p bed {}; fi' ::: /g/data/<HPC_PROJECT>/data/atac/atac_count_outs/force_cells_test/*/*/outs/fragments_nostrand.tsv.gz
