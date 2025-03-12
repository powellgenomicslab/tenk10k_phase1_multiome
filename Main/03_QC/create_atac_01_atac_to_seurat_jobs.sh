####
# folders=($(find /g/data/fy54/data/atac/atac_count_outs/ -type d -name "S0*" -exec basename {} \; | sort -V))

for i in {220..221}
do

sed 's/pool/'S0${i}'/g' atac_01_atac_to_seurat_pool_1.qsub.sh > atac_01_atac_to_seurat_S0${i}_1.qsub.sh
sed 's/pool/'S0${i}'/g' atac_01_atac_to_seurat_pool_2.qsub.sh > atac_01_atac_to_seurat_S0${i}_2.qsub.sh
qsub atac_01_atac_to_seurat_S0${i}_1.qsub.sh
qsub atac_01_atac_to_seurat_S0${i}_2.qsub.sh
rm atac_01_atac_to_seurat_S0${i}_1.qsub.sh
rm atac_01_atac_to_seurat_S0${i}_2.qsub.sh

done



####
