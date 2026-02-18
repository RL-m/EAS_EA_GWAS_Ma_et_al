refFile=$1
sumFile=$2
outFile=$3
Plink19="/plink-1.9/plink"
${Plink19} \
--bfile ${refFile} \
--clump ${sumFile} \
--clump-p1 5e-8 \
--clump-p2 5e-8 \
--clump-r2 0.1 \
--clump-kb 10000 \
--out ${outFile}