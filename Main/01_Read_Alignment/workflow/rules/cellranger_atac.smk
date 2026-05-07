# ⚙️ functions 

def get_target_cell_count(wildcards):
    """
    Get the target number of cells for the capture.
    This will be used in the cellranger count --force-cells parameter
    """
    
    # read in the lab metadata (see workflow/scripts/01_generate_sample_sheets.py)
    master_df = pd.read_csv('/g/data/<HPC_PROJECT>/repos/tenk10k_cellranger/sample_sheets/mastertable.csv', index_col=0)
    sample_target_cells = master_df.loc[master_df["sample_id_cellranger"] == wildcards.sample, ["Targeted Cell Number"]]

    if sample_target_cells.empty:
        raise ValueError(f"Sample '{wildcards.sample}' not found in the metadata.")
    if len(sample_target_cells) > 1:
        raise ValueError(f"Multiple entries found for sample '{wildcards.sample}'. Expected only one.")
    
    return sample_target_cells.squeeze()

# 📏 rules

rule cellranger_atac_count:
    input:
        fastqs=f"{config['fastqs']}/{{sequencing_date}}_{{flowcell_id}}/{{flowcell_id}}",
        reference=config["reference"]
    output:
        outs_dir=directory(os.path.join(config["atac_count_outdir"], "{sequencing_date}_{flowcell_id}", "{sample}/outs")),
        command=os.path.join(config["atac_count_outdir"], "{sequencing_date}_{flowcell_id}", "{sample}/outs/_command"),
        web_summary=os.path.join(config["atac_count_outdir"], "{sequencing_date}_{flowcell_id}", "{sample}/outs/web_summary.html"),
        web_summary_copy=os.path.join(config["atac_count_outdir"], "web_summaries", "{sample}_{sequencing_date}_{flowcell_id}_web_summary.html")
    params:
        cellranger_atac=config["cellranger_atac"],
        atac_count_outdir_flowcell=os.path.join(config["atac_count_outdir"],"{sequencing_date}_{flowcell_id}"),
        target_cells=get_target_cell_count
    conda:
        "cellranger-atac-snakemake"
    resources:
        localcores=48, # each job uses max resources for a for cascade lake node 
        localmem=190,
        queue="normal",
        job_storage_gb=400, # maxxing out the node any way so I don't think its any cheaper to specify less storage 
        walltime="20:00:00",
        iointensive=2 # number of i/o intensive volumes storage (1 TB storage each) 
    log: 
        count_log="logs/cellranger_atac_count/count_{sample}_{sequencing_date}_{flowcell_id}_snakemake.log"
    shell:
        """
        (
            # NOTE: Snakemake will automatically create the output sample directory HOWEVER, cellranger will error unless the output directory is created by cellranger
            #     my solution is to run in tmpdir and then move the output files across to the pipeline output dir.
            #     running in tempdir allows less storage usage by the large intermediate files

            cd /iointensive
  
            CMD="{params.cellranger_atac} count \
                --id={wildcards.sample} \
                --sample={wildcards.sample} \
                --reference={input.reference} \
                --fastqs={input.fastqs}/{wildcards.sample} \
                --localcores={resources.localcores} \
                --localmem={resources.localmem} \
                --force-cells={params.target_cells} \
                --disable-ui"
            
            echo $CMD 
            
            eval $CMD 

            echo $CMD > {wildcards.sample}/outs/_command
            
            rm -rf {wildcards.sample}/SC_ATAC_COUNTER_CS/

            # make output directory 
            if [ ! -d {params.atac_count_outdir_flowcell} ]; then
                mkdir -p {params.atac_count_outdir_flowcell}
            fi

            # move output from tempdir into target dir, overwrite if previous output exists 
            rm -rf {params.atac_count_outdir_flowcell}/{wildcards.sample} 
            mv {wildcards.sample} {params.atac_count_outdir_flowcell}

            cp {output.web_summary} {output.web_summary_copy}

        ) &> {log.count_log}

        """