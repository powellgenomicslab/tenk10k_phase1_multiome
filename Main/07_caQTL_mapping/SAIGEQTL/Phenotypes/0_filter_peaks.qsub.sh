## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -o /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Logs/0_filter_peaks.OU
#$ -e /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes/Logs/0_filter_peaks.ER
#$ -N filter_peaks
#$ -q short.q
#$ -pe smp 16
#$ -l mem_requested=30G
#$ -l tmp_requested=200G
#$ -r yes

cd /directflow/SCCGGroupShare/projects/oscdon/SAIGEQTL/Phenotypes
conda activate py

touch Logs/0_filter_peaks.log
> Logs/0_filter_peaks.log

# Main command and record the time and memory usage
python -u 0_filter_peaks.py >> Logs/0_filter_peaks.log 2>&1

####