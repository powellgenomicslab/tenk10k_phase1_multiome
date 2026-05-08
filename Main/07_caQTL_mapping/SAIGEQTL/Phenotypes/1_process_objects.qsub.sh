## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -o /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Logs/1_process_objects.OU
#$ -e /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Logs/1_process_objects.ER
#$ -N filter_peaks
#$ -q short.q
#$ -pe smp 1
#$ -l mem_requested=32G
#$ -l tmp_requested=20G
#$ -r yes
#$ -t 1-231

# Extract pool name
i=${SGE_TASK_ID};
POOLS=/directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/OldRun/pool_list.txt
POOL=$(sed -n "${i}p" "$POOLS")

cd /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes
conda activate py

touch Logs/1_process_objects/${POOL}.log
> Logs/1_process_objects/${POOL}.log

# Main command and record the time and memory usage
python -u 1_process_objects.py ${POOL} >> Logs/1_process_objects/${POOL}.log 2>&1

####