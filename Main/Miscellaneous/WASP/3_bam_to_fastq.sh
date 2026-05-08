#!/bin/bash
#PBS -P fy54
#PBS -q normalbw
#PBS -N 3_bam_to_fastq
#PBS -l walltime=4:00:00
#PBS -l storage=gdata/ei56+gdata/fy54
#PBS -l mem=512GB
#PBS -l jobfs=100GB
#PBS -l ncpus=112
#PBS -r y
#PBS -l wd
#PBS -o /g/data/fy54/od8037/TenK10K/WASP/OutputLogs/3_bam_to_fastq/3_bam_to_fastq.OU
#PBS -e /g/data/fy54/od8037/TenK10K/WASP/OutputLogs/3_bam_to_fastq/3_bam_to_fastq.ER

bam_to_fastq() {
    local POOL=${1}
    local CELLTYPE=${2}
    local DONOR=${3}

    local OUTDIR=/g/data/fy54/od8037/TenK10K/WASP/FASTQs/${POOL}-${CELLTYPE}-${DONOR}
    if [ -d "${OUTDIR}" ];then
        rm -rf "${OUTDIR}"
    fi

    /g/data/ei56/od8037/software/cellranger-atac-2.2.0/cellranger-atac bamtofastq \
        --relaxed \
        --traceback \
        --nthreads=2 \
        /g/data/fy54/od8037/TenK10K/WASP/BAMs/${POOL}-${CELLTYPE}-${DONOR}.bam \
        ${OUTDIR}
} 
export -f bam_to_fastq

module load parallel
parallel -j 56 --colsep '\t' --verbose \
    bam_to_fastq {1} {2} {3} \
    :::: /g/data/fy54/od8037/TenK10K/WASP/combinations.tsv

echo "Finished"
