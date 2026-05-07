#!/bin/bash
#PBS -P <HPC_PROJECT>
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/<HPC_PROJECT>+gdata/<HPC_PROJECT_2>
#PBS -l mem=16GB
#PBS -l ncpus=1
#PBS -l jobfs=16GB
#PBS -l wd
#PBS -e /g/data/<HPC_PROJECT>/repos/tenk10k_cellranger/logs/pbs_logs/snakemake_runner.stderr
#PBS -o /g/data/<HPC_PROJECT>/repos/tenk10k_cellranger/logs/pbs_logs/snakemake_runner.stdout
#PBS -M <YOUR_EMAIL>
#PBS -m ae

# USAGE: Before running, create a config.yaml and specify the input bcl dir and sample sheet in the config.
#        Update $CONFIG to point to the config.yaml for the run you want to execute.
#        This version of the runner uses --force-cells which is determined by the targeted number of cells specified in the library metadata on airtable 
#        This is intended to be run on one flowcell at a time

# Example:

CONFIG=/g/data/<HPC_PROJECT>/repos/tenk10k_cellranger/config_files/config_files_force_cells/<DATE>_<INSTRUMENT>_<RUN_NUMBER>_<FLOWCELL_ID>.yaml

# Do not change the following:

SNAKEMAKE=<PATH_TO_SNAKEMAKE>
SNAKEFILE="${PBS_O_WORKDIR}/Snakefile"

RUN_NAME=$(basename ${CONFIG%_config.yaml})

(
    ${SNAKEMAKE} --snakefile ${SNAKEFILE} \
        --configfile ${CONFIG} \
        --jobs 100 \
        --use-conda \
        --executor cluster-generic \
        --cluster-generic-submit-cmd \
            "qsub -S /bin/bash \
                -N {rule} \
                -q {resources.queue} \
                -P <HPC_PROJECT> \
                -l mem={resources.localmem}GB \
                -l ncpus={resources.localcores} \
                -l jobfs={resources.job_storage_gb}GB \
                -l walltime={resources.walltime} \
                -l storage=gdata/<HPC_PROJECT>+gdata/<HPC_PROJECT_2>+scratch/<HPC_PROJECT> \
                -l iointensive={resources.iointensive} \
                -l wd \
                -M <YOUR_EMAIL> \
                -m ae \
                -e /g/data/<HPC_PROJECT>/repos/tenk10k_cellranger/logs/pbs_logs/snakemake_{rule}_{wildcards.sequencing_date}.stderr \
                -o /g/data/<HPC_PROJECT>/repos/tenk10k_cellranger/logs/pbs_logs/snakemake_{rule}_{wildcards.sequencing_date}.stdout" --rerun-incomplete

) &> "/g/data/<HPC_PROJECT>/repos/tenk10k_cellranger/logs/pipeline/${RUN_NAME}.out"


# -n for dry run
# -F for force rerun all
