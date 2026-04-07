#!/bin/bash
# Request SGE resources
#$ -N BH_DEM4            # Name of the job
#$ -cwd                 # Run job from the current directory
#$ -l mem_requested=8G # Memory requested per core
#$ -pe smp 24            # Request 4 CPU cores (adjust based on your needs)
#$ -l h_vmem=8G        # Set hard memory limit per core
#$ -l tmp_requested=8G # Request 20 GB of temporary space
#$ -r yes
#$ -o output4.log        # Log file for standard output
#$ -e error4.log         # Log file for error output

source ~/miniconda3/etc/profile.d/conda.sh

# Load R module (if your HPC uses modules)
conda activate scenicplus

SIF=/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/BioHeart_LBIO/resource/Demuxafy.sif
BIND_PATH=/directflow/SCCGGroupShare/projects

THREADS=24
N=16

INDEX=ES001_P4_4
NOCHR_VCF=/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/BioHeart_LBIO/output/library_extract/BioHeart/${INDEX}.vcf.gz
VCF=/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/BioHeart_LBIO/output/library_extract/BioHeart/${INDEX}_chrmap.vcf.gz
VCFTEM=/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/BioHeart_LBIO/output/library_extract/BioHeart/${INDEX}_chrmap_temp.vcf.gz
# DONOR_ID=/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/BioHeart_LBIO/output/library_extract/BioHeart/ES001_P1_1_samplelist.txt

BARCODES=/directflow/SCCGGroupShare/projects/josealq/multiome/data/ATAC/${INDEX}/outs/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
BAM=/directflow/SCCGGroupShare/projects/josealq/multiome/data/ATAC/${INDEX}/outs/outs/atac_possorted_bam.bam

FASTA=/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/BioHeart_LBIO/resource/genome.fa
VIREO_OUTDIR=/directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/BioHeart_LBIO/output/demultiplexing/BioHeart/${INDEX}
mkdir -p $VIREO_OUTDIR

# Transform chromatin label from '1', '2', '3' -> 'chr1', 'chr2', 'chr3'
echo "Step 0"
bcftools annotate --rename-chrs /directflow/SCCGGroupShare/projects/jayfan/Projects/Multiome/tenk10k_phase1/BioHeart_LBIO/resource/chr_map.txt $NOCHR_VCF > $VCF

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
singularity exec --bind $BIND_PATH $SIF bgzip -c $VCF > $VCFTEM
rm $VCF
mv $VCFTEM $VCF
singularity exec --bind $BIND_PATH $SIF tabix -p vcf $VCF

# Match SNP and donor list
echo "Step 3"
singularity exec --bind $BIND_PATH $SIF bcftools view $VCF \
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

