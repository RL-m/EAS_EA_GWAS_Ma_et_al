chr=$1
bFile=$2
maFile=$3
outFile=$4
gctaAdd="/packages/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
geneList="glist_ensgid_hg38_v40.txt"
# Output mBAT-combo, mBAT and fastBAT p-values:
${gctaAdd} --bfile ${bFile} \
--mBAT-combo ${maFile} \
--mBAT-gene-list ${geneList} \
--mBAT-print-all-p \
--chr ${chr} \
--out ${outFile} \
--thread-num 10