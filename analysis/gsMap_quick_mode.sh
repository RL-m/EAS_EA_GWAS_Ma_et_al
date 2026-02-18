mainAdd=$1
sumstats=$2
trait=$3
gsmap quick_mode \
    --workdir ${mainAdd} \
    --homolog_file '/gsMap/gsMap_resource/homologs/mouse_human_homologs.txt' \
    --sample_name 'Adult_Mouse_Brain' \
    --gsMap_resource_dir '/gsMap/gsMap_resource' \
    --hdf5_path '/gsMAP/Adult_Mouse_brain.h5ad' \
    --annotation 'annotation' \
    --data_layer 'count' \
    --sumstats_file ${sumstats} \
    --trait_name ${trait}