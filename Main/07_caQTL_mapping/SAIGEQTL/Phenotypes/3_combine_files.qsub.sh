## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -o /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Logs/3_combine_files/celltype.OU
#$ -e /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Logs/3_combine_files/celltype.ER
#$ -N celltype
#$ -q short.q
#$ -pe smp 8
#$ -l mem_requested=40G
#$ -l tmp_requested=30G
#$ -r yes
#$ -t 1-22

# Extract chromosome number
CHR_NUM=${SGE_TASK_ID};

cd /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes
conda activate py

touch Logs/3_combine_files/celltype/chr${CHR_NUM}.log
> Logs/3_combine_files/celltype/chr${CHR_NUM}.log

# Main command and record the time and memory usage
python -u 3_combine_files.py celltype ${CHR_NUM} >> Logs/3_combine_files/celltype/chr${CHR_NUM}.log 2>&1

####