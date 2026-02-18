# snpLoc: SNP,CHR,POS #
# snpP: SNP,P #
snpLoc=$1
snpP=$2
outFile=$3
N=$4
refFile=$5
magma="/packages/MAGMA/magma"
geneLoc="/packages/MAGMA/NCBI38.gene.loc"
${magma} \
--annotate \
--snp-loc ${snpLoc} --gene-loc ${geneLoc} \
--out ${outFile}
### gene analysis ###
${magma} \
--bfile ${refFile} \
--pval ${snpP} N=${N} \
--gene-annot ${outFile}.genes.annot \
--out ${outFile}
### gene set analysis ###
setAnnot="/packages/MAGMA/msigdb.v2023.2.Hs.entrez.gmt.txt"
${magma} \
--gene-results ${outFile}.genes.raw \
--set-annot ${setAnnot} \
--out ${outFile}