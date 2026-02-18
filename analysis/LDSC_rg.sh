sumFile1=$1
sumFile2=$2
outFile=$3
source ~/.bashrc
source /packages/miniconda3/etc/profile.d/conda.sh
conda activate ldsc
ldsc="/packages/ldsc/ldsc.py"
weight="/LDSC_resource/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC."
${ldsc} \
--rg ${sumFile1},${sumFile2} \
--ref-ld-chr ${weight} \
--w-ld-chr ${weight} \
--out ${outFile}