source ~/.bashrc
source-conda
conda activate ldsc
inFile=$1
outFile=$2
ldsc="/packages/ldsc/ldsc.py"
weight="/LDSC_resource/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC."
${ldsc} \
--h2 ${inFile} \
--ref-ld-chr ${weight} \
--w-ld-chr ${weight} \
--out ${outFile}