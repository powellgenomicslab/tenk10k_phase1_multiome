
import anndata
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import os
import scglue

# %%
rna = anndata.read_h5ad("/g/data/fy54/jf1058/ei56_nomem/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/10XMultiome/Revision/10Xdata_cell_type_spec_CRISPR/CD14_Mono_enriched_gene/save/s01_preprocessing/omics_data/rna.h5ad")
atac = anndata.read_h5ad("/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/data_processing/cell_type_spec_TenK10K_phase1/ATAC_classification/output_celltype_238lib/CD14_Mono/ATAC_CD14_Mono.h5ad")

split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)

rna.var["name"] = rna.var_names
atac.var["name"] = atac.var_names

# %%
genes = scglue.genomics.Bed(rna.var.assign(name=rna.var_names))
peaks = scglue.genomics.Bed(atac.var.assign(name=atac.var_names))
tss = genes.strand_specific_start_site()
promoters = tss.expand(2000, 0)
flanks = tss.expand(500, 500)

# %%
def aggregate_SMR(dir):
    for chr in range(1, 23):
        if os.path.exists(dir + "Chr" + str(chr) + "_results.smr") == False:
            continue
        SMR_chr = pd.read_csv(dir + "Chr" + str(chr) + "_results.smr", sep = "\t")
        if chr == 1:
            SMR_agg = SMR_chr
        else:
            SMR_agg = pd.concat([SMR_agg, SMR_chr])
    return SMR_agg
smr_peak_gene = aggregate_SMR("/g/data/fy54/jf1058/ei56_nomem/tenk10k_phase1/SMR/Revision_flipping_corrected_genotypes/output/caQTL2eQTL_1Mb/CD14_Mono/")
smr_sig_peak_gene = smr_peak_gene[smr_peak_gene['p_SMR'] < 5e-8]

chip = scglue.genomics.read_bed("/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/10XMultiome/data/ENCODE-TF-ChIP-hg38.bed.gz")

######################################

GLUE_scores = pd.read_pickle("/g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/TenK10K_231lib_celltype_eQTL/CD14_Mono/save/s04_infer_gene_tf/gene_peak_conn.pkl.gz")
GLUE_scores = GLUE_scores[GLUE_scores["glue"] > 0.8]

gene_of_interest = "APOBEC3A"
GLUE_scores_GOI = GLUE_scores[GLUE_scores['gene'] == gene_of_interest][['gene', 'peak']]
smr_scores_GOI = smr_sig_peak_gene[smr_sig_peak_gene['Outco_Gene'] == gene_of_interest][['Outco_Gene', 'Expo_ID']].rename(columns={'Outco_Gene': 'gene', 'Expo_ID': 'peak'})

SAVE_PATH = "/g/data/fy54/jf1058/ei56_nomem/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/TenK10K_231lib_celltype_eQTL/CD14_Mono/save/s21_pyGenomeTrack_GLUE_SMR_APOBEC3A"
os.makedirs(SAVE_PATH, exist_ok = True)

scglue.genomics.Bed(atac.var).write_bed(f"{SAVE_PATH}/peaks.bed", ncols=3)

gene2peak_link_smr = smr_scores_GOI.merge(
        scglue.genomics.Bed(rna.var).strand_specific_start_site().df.iloc[:, :4], how="left", left_on="gene", right_on="name"
    ).merge(
        scglue.genomics.Bed(atac.var).df.iloc[:, :4], how="left", left_on="peak", right_on="name"
    ).loc[:, [
        "chrom_x", "chromStart_x", "chromEnd_x",
        "chrom_y", "chromStart_y", "chromEnd_y"
    ]].dropna().assign(score=1)

gene2peak_link_smr["chromStart_x"] = gene2peak_link_smr["chromStart_x"].astype('int')
gene2peak_link_smr["chromEnd_x"] = gene2peak_link_smr["chromEnd_x"].astype('int')
gene2peak_link_smr["chromStart_y"] = gene2peak_link_smr["chromStart_y"].astype('int')
gene2peak_link_smr["chromEnd_y"] = gene2peak_link_smr["chromEnd_y"].astype('int')
gene2peak_link_smr.to_csv(f"{SAVE_PATH}/gene2peak_smr_{gene_of_interest}.links", sep="\t", index=False, header=False)

# %%
# Refine the GLUE score, only focus on the SMR highlighted gene-peak pairs
gene2peak_link_glue = GLUE_scores_GOI.merge(
       scglue.genomics.Bed(rna.var).strand_specific_start_site().df.iloc[:, :4], how="left", left_on="gene", right_on="name"
   ).merge(
       scglue.genomics.Bed(atac.var).df.iloc[:, :4], how="left", left_on="peak", right_on="name"
   ).loc[:, [
       "chrom_x", "chromStart_x", "chromEnd_x",
       "chrom_y", "chromStart_y", "chromEnd_y"
   ]].dropna()
gene2peak_link_glue["chromStart_x"] = gene2peak_link_glue["chromStart_x"].astype('int')
gene2peak_link_glue["chromEnd_x"] = gene2peak_link_glue["chromEnd_x"].astype('int')
gene2peak_link_glue["chromStart_y"] = gene2peak_link_glue["chromStart_y"].astype('int')
gene2peak_link_glue["chromEnd_y"] = gene2peak_link_glue["chromEnd_y"].astype('int')
gene2peak_link_glue.to_csv(f"{SAVE_PATH}/gene2peak_glue_{gene_of_interest}.links", sep="\t", index=False, header=False)

# Create another link file for ATAC peak distance
tss_chrom, tss_pos = 'chr22', 38952740
window = 10_000  # +/- 10 kb

same_chr = pd.DataFrame(peaks[peaks['chrom'] == tss_chrom]).copy()
same_chr['distance'] = ((same_chr['chromStart'] - tss_pos).clip(lower=0)
                      + (tss_pos - same_chr['chromEnd']).clip(lower=0))

nearby = same_chr[same_chr['distance'] <= window].reset_index(drop=True)

result = pd.DataFrame({
    'chrom_x':      tss_chrom,
    'chromStart_x': tss_pos,
    'chromEnd_x':   tss_pos + 1,
    'chrom_y':      nearby['chrom'].values,
    'chromStart_y': nearby['chromStart'].values,
    'chromEnd_y':   nearby['chromEnd'].values,
})
result.to_csv(f"{SAVE_PATH}/gene2peak_dist_{gene_of_interest}.links",
              sep="\t", index=False, header=False)


# Modify tracks_SMR.ini file here, change the glue score links file.
subprocess.run('pyGenomeTracks --tracks tracks_revision_APOBEC3A.ini --region chr22:38940000-39000000 --outFileName /g/data/fy54/jf1058/ei56_nomem/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/TenK10K_231lib_celltype_eQTL/CD14_Mono/save/s21_pyGenomeTrack_GLUE_SMR_APOBEC3A/tracks_revision_APOBEC3A.png --dpi 300', shell = True, executable="/bin/bash")

zcat /g/data/ei56/jf1058/TEMP/Brenner_copy/tenk10k_phase1/scGLUE/metadata/gencode.v47.chr_patch_hapl_scaff.annotation.gtf.gz \
  | awk '/^#/ || /gene_type "protein_coding"/' \
  | gzip > /g/data/fy54/jf1058/ei56_nomem/tenk10k_phase1/scGLUE/src_ATAC_Manuscript/TenK10K_231lib_celltype_eQTL/CD14_Mono/save/s13_pyGenomeTrack_GLUE_SMR/genes.protein_coding.gtf.gz

