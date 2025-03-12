#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N vireo_demultiplex_TOB_ATAC_pool_1
#PBS -l walltime=5:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=96GB
#PBS -l ncpus=24
#PBS -r y
#PBS -l wd
#PBS -M a.xue@garvan.org.au
#PBS -m ae

DIR=/g/data/ei56/ax3061/proj/tenk10k/caQTL/demultiplexing/output/Vireo
SIF=/g/data/ei56/ax3061/software/Demuxafy.sif
BIND_PATH=/g/data/ei56/ax3061

THREADS=24
N=8

VCF=/g/data/ei56/projects/TenK10K/genotype/pool.vcf
DONOR_ID=/g/data/ei56/projects/TenK10K/genotype/pool_sample_list.txt
DATA_DIR=/g/data/fy54/data/atac/atac_count_outs
folders=($(find ${DATA_DIR} -maxdepth 2 -type d -name "pool_1*" -exec readlink -f {} \;))
# Check if the folders array is empty
if [ ${#folders[@]} -eq 0 ]; then
    folders=($(find "${DATA_DIR}" -maxdepth 2 -type d -name "pool_R_1" -exec readlink -f {} \;))
fi
echo ${folders}
BARCODES=${folders}/outs/filtered_peak_bc_matrix/barcodes.tsv
BAM=${folders}/outs/possorted_bam.bam

FASTA=/g/data/ei56/projects/TenK10K/reference/genome.fa
VIREO_OUTDIR=/g/data/ei56/ax3061/proj/tenk10k/caQTL/demultiplexing/output/Vireo/pool_1
mkdir -p $VIREO_OUTDIR

module load singularity

# CellSNP Pileup
echo "Step 1"
singularity exec --bind $BIND_PATH $SIF cellsnp_pileup.py \
          -s $BAM \
          -b $BARCODES \
          -O $VIREO_OUTDIR \
          -R $VCF \
          --cellTAG CB \
          --UMItag None \
          -p $THREADS \
          --minMAF 0.1 \
          --minCOUNT 20 \
          --gzip GZIP

# bgzip and index the files
echo "Step 2"
singularity exec --bind $BIND_PATH $SIF bgzip -c $VCF > $VCF.gz
singularity exec --bind $BIND_PATH $SIF tabix -p vcf $VCF.gz

# Match SNP and donor list
echo "Step 3"
singularity exec --bind $BIND_PATH $SIF bcftools view $VCF.gz \
        -R $VIREO_OUTDIR/cellSNP.base.vcf.gz \
        -Ov \
        -o $VIREO_OUTDIR/donor_subset.vcf

# Run Vireo
echo "Step 4"
singularity exec --bind $BIND_PATH $SIF vireo \
        -c $VIREO_OUTDIR \
        -d $VIREO_OUTDIR/donor_subset.vcf \
        -o $VIREO_OUTDIR \
        -t GT \
        -N $N \
        -p $THREADS \
        --callAmbientRNAs



