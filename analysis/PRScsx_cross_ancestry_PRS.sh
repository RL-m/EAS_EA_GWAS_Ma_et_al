# Prepare the summary file in the name of ${pop1}.txt and ${pop2}.txt in the ${outAdd} folder #
# sumFile: SNP,A1,A2,BETA,SE #
chr=${TASK_ID}
pop1=EAS
pop2=EUR
N1=$1
N2=$2
refFile="/PRScsx_reference"
targetGeno=$3
outAdd=$4
PRScsx="/packages/PRScsx/PRScsx.py"
python ${PRScsx} \
--ref_dir=${refFile} --bim_prefix=${targetGeno} \
--sst_file=${outAdd}/${pop1}.txt,${outAdd}/${pop2}.txt --meta=TRUE \
--n_gwas=${N1},${N2} --pop=${pop1},${pop2} --chrom=${chr} \
--out_dir=${outAdd} --out_name=PRScsx.${pop1}.${pop2}.txt