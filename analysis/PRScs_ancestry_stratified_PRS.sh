# run PRScs #
# sumFile: SNP,A1,A2,BETA,SE #
chr=${TASK_ID}
refFile="/PRScsx_reference/ldblk_1kg_eas"
# refFile="/PRScsx_reference/ldblk_1kg_eur"
targetGeno=$1
sumFile=$2
N=$3
outFile=$4
PRScs=/packages/PRScs/PRScs.py
python ${PRScs} \
--ref_dir=${refFile} \
--bim_prefix=${targetGeno} \
--sst_file=${sumFile} \
--n_gwas=${N} \
--chrom=${chr} \
--out_dir=${outFile}

# construct PRS by PLINK #
Plink2Add="/packages/plink2"
targetGeno=$1
snpEffect=$2
prsFile=$3
${Plink2Add} \
--bfile ${targetGeno} \
--score ${snpEffect} 2 4 6 header \
--out ${prsFile}