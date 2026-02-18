bFile=$1
exposureList=$2
outcomeList=$3
outFile=$4
GCTA="/packages/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
${GCTA} \
--bfile ${bFile} \
--gsmr-file ${exposureList} ${outcomeList} \
--gsmr-direction 2 --effect-plot --thread-num 4 \
--out ${outFile}