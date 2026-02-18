sumstats=$1
outFile=$2
refAdd="/LDSC_resource"
baseline=${refAdd}/EAS_baselineLD_v2.2/baselineLD.
weight=${refAdd}/1000G_Phase3_EAS_weights_hm3_no_MHC/weights.EAS.hm3_noMHC.
frq=${refAdd}/1000G_Phase3_frq/1000G.EAS.QC.
source ~/.bashrc
source-conda
conda activate ldsc
ldsc="/packages/ldsc/ldsc.py"
python ${ldsc} \
--h2 ${sumstats} \
--ref-ld-chr ${baseline} \
--w-ld-chr ${weight} \
--overlap-annot \
--frqfile-chr ${frq} \
--print-coefficients \
--out ${outFile}