
# ⚙️ functions 

# 📏 rules

rule cellranger_mkfastq:
    input:
        bcl_dir=config["bcl_dir"],
        samplesheet=os.path.join(config["sample_sheet_dir"], "sample_sheet_{sequencing_date}_{flowcell_id}.csv")
    output:
        fastqs=directory(f"{config['fastqs']}/{{sequencing_date}}_{{flowcell_id}}/{{flowcell_id}}"),
        out_dir=directory(f"{config['fastqs']}/{{sequencing_date}}_{{flowcell_id}}"),
        command=f"{config['fastqs']}/{{sequencing_date}}_{{flowcell_id}}/_command"
    params:
        cellranger_atac=config["cellranger_atac"],
    conda:
        "bcl2fastq-env"
    log:
        mkfastq_log="logs/mkfastq/mkfastq_snakemake_{sequencing_date}_{flowcell_id}.log"
    resources:
        localcores=48,
        localmem=190,
        queue="normal",
        job_storage_gb=400,
        walltime="20:00:00",
        iointensive=0

    shell:
        """
        (
            cd ${{TMPDIR}}
            mkdir -p {output.out_dir}

            CMD="{params.cellranger_atac} mkfastq \
                --id=mkfastq \
                --run={input.bcl_dir} \
                --csv={input.samplesheet} \
                --output-dir={output.out_dir} \
                --jobmode=local \
                --localcores={resources.localcores} \
                --localmem={resources.localmem} \
                --disable-ui"
            echo $CMD > {output.command}
            eval $CMD

        ) &> {log.mkfastq_log}
        """
