sumFile=$1
refFile=$2
qtlFile=$3
outFile=$4
smrAdd="/packages/smr_Linux"
${smrAdd} \
--bfile ${refFile} \
--gwas-summary ${sumFile} \
--beqtl-summary ${qtlFile} \
--maf 0.01 --smr-multi --thread-num 2 \
--out ${outFile}