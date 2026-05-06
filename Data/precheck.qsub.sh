## SGE SETTINGS
#$ -cwd
#$ -S /bin/bash
#$ -o stdout_precheck
#$ -e stderr_precheck
#$ -N precheck
#$ -q short.q
#$ -pe smp 1
#$ -l mem_requested=60G
#$ -r yes


cd /directflow/SCCGGroupShare/projects/angxue/proj/multiome/

Rscript421=/share/ScratchGeneral/angxue/software/miniconda3/envs/R421/bin/Rscript

conda activate R421

touch precheck.log
$Rscript421 precheck.R >> precheck.log


####
