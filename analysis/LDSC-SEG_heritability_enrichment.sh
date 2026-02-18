cts_name=$1
sumstats=$2
outFile=$3
source ~/.bashrc
source-conda
conda activate ldsc
ldsc="/packages/ldsc/ldsc.py"
refAdd="/LDSC_resource"
baseline=${refAdd}/EAS_baselineLD_v2.2/baselineLD.
weight=${refAdd}/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC.
ldcts=${refAdd}/LDSC-SEG/${cts_name}.ldcts
${ldsc} \
--h2-cts ${sumstats} \
--ref-ld-chr ${baseline} \
--out ${outFile} \
--ref-ld-chr-cts ${ldcts} \
--w-ld-chr ${weight}