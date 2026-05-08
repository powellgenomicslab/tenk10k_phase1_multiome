#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N 4_Part_num
#PBS -l walltime=24:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=128GB
#PBS -l jobfs=100GB
#PBS -l ncpus=28
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/WASP/OutputLogs/4_star_wasp/Part_num.OU
#PBS -e /g/data/fy54/od8037/TenK10K/WASP/OutputLogs/4_star_wasp/Part_num.ER

export TMPDIR
run_star_wasp() {
    local POOL=${1}
    local CELLTYPE=${2}
    local DONOR=${3}

    local OUTDIR=${TMPDIR}/${CELLTYPE}_${DONOR}
    local FINALDIR=/g/data/fy54/od8037/TenK10K/WASP/Outputs/FilteredBAMs

    # Run STAR with WASP mapping bias correction
    cd ${TMPDIR}

    /g/data/fy54/od8037/software/STAR-2.7.11b/STAR/bin/Linux_x86_64_static/STAR \
        --runThreadN 7 \
        --outFileNamePrefix ${OUTDIR}/out_ \
        --outTmpDir ${TMPDIR}/TMP_${CELLTYPE}_${DONOR} \
        --genomeDir /g/data/fy54/od8037/TenK10K/WASP/Genome \
        --soloCBwhitelist /g/data/fy54/od8037/TenK10K/WASP/FromBrenner/737K-cratac-v1.txt \
        --readFilesPrefix /g/data/fy54/od8037/TenK10K/WASP/FASTQs/${POOL}-${CELLTYPE}-${DONOR}/*/ \
        --readFilesIn bamtofastq_S1_L001_R1_001.fastq.gz,bamtofastq_S1_L002_R1_001.fastq.gz bamtofastq_S1_L001_R3_001.fastq.gz,bamtofastq_S1_L002_R3_001.fastq.gz bamtofastq_S1_L001_R2_001.fastq.gz,bamtofastq_S1_L002_R2_001.fastq.gz \
        --readFilesCommand zcat \
        --readFilesSAMattrKeep None \
        --soloType CB_samTagOut \
        --soloCBmatchWLtype 1MM \
        --soloBarcodeReadLength 0 \
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNmax 4 \
        --alignIntronMax 1 \
        --alignMatesGapMax 1000 \
        --outSAMtype BAM Unsorted \
        --outSAMattributes NH HI AS nM CB sS sQ ha \
        --waspOutputMode SAMtag \
        --varVCFfile /g/data/fy54/od8037/TenK10K/WASP/VCFs/DonorVCFs/${DONOR}.vcf

    # Filter out reads that did not pass WASP filtering (reads with vW:i:1 are passes, all other values for vW are fails)
    /g/data/ei56/od8037/software/samtools/bin/samtools view \
        -e '[vW] == 1 || !exists([vW])' \
        -o ${OUTDIR}/Aligned.out.filtered.bam ${OUTDIR}/out_Aligned.out.bam \
        --threads 7

    # Sort and index the filtered BAM
    /g/data/ei56/od8037/software/samtools/bin/samtools sort ${OUTDIR}/Aligned.out.filtered.bam \
        -o ${FINALDIR}/${CELLTYPE}-${DONOR}.bam \
        --threads 7 
    /g/data/ei56/od8037/software/samtools/bin/samtools index ${FINALDIR}/${CELLTYPE}-${DONOR}.bam  

    rm -rf ${TMPDIR}/${CELLTYPE}_${DONOR}
    rm -rf ${TMPDIR}/TMP_${CELLTYPE}_${DONOR}
}
export -f run_star_wasp

module load parallel
parallel -j 4 --colsep '\t' --verbose \
    run_star_wasp {1} {2} {3} \
    :::: /g/data/fy54/od8037/TenK10K/WASP/Combinations/combinations_num.tsv

echo "Finished"
