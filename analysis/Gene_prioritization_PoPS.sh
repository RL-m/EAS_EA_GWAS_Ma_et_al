source ~/.bashrc
source-conda
conda activate pops
PoPS="/packages/pops-master/pops.py"
gene_annot="/packages/pops-master/example/data/utils/gene_annot_jun10.txt"
featureFile="/packages/pops-master/pops_all_features_munged/all_pops_features"
control_feature="/packages/pops-master/example/data/utils/features_jul17_control.txt"
filePrefix=$1
outFile=${filePrefix}.PoPS
magmaResult=${filePrefix}.MAGMA
python ${PoPS} \
--verbose \
--gene_annot_path ${gene_annot} \
--feature_mat_prefix ${featureFile}  \
--num_feature_chunks 116 \
--magma_prefix ${magmaResult} \
--control_features_path ${control_feature} \
--out_prefix ${outFile}