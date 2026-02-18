refFile=$1
maFile=$2
outFile=$3
gctaAdd="/packages/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
${gctaAdd} \
--bfile ${refFile} \
--cojo-p 5e-08 \
--cojo-file ${maFile} \
--cojo-slct --threads 2 \
--out ${outFile}