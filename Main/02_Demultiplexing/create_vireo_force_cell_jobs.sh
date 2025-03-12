####
# folders=($(find /g/data/fy54/data/atac/atac_count_outs/ -type d -name "S0*" -exec basename {} \; | sort -V))

for i in {304..338}
do
sed 's/pool/'S0${i}'/g' force_cell_vireo_TOB_ATAC_pool_1.qsub.sh > force_cell_vireo_TOB_ATAC_S0${i}_1.qsub.sh
sed 's/pool/'S0${i}'/g' force_cell_vireo_TOB_ATAC_pool_2.qsub.sh > force_cell_vireo_TOB_ATAC_S0${i}_2.qsub.sh
qsub force_cell_vireo_TOB_ATAC_S0${i}_1.qsub.sh
qsub force_cell_vireo_TOB_ATAC_S0${i}_2.qsub.sh
done



####
