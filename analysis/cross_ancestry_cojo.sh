refFile="/1000G/GRCh38_EAS/1kGP_high_coverage_Illumina.EAS_unrelated.geno05.mind05.maf01.hwe1e6.autosome"
sumFile="/Meta/EAS_EA_GWAS.ma"
snpList="EA4_reported_clumping_SNP.extension_by_UKB50k_ld0.1.addMissing.snplist"
outFile="EA4_sigSNPs_in_EAS_EA_GWAS"
gctaAdd="/packages/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
${gctaAdd}  \
--bfile ${refFile} \
--cojo-file ${sumFile} \
--extract ${snpList} \
--cojo-slct \
--threads 4 \
--out ${outFile}

sumFile="/Meta/EAS_EA_GWAS.ma"
snpList="EA4_sigSNPs_in_EAS_EA_GWAS.snplist"
sed '1d' /EA4_sigSNPs_in_EAS_EA_GWAS.jma.cojo | cut -f2 > ${snpList}
outFile="EAS_EA_GWAS.after_cond_EA4_sigSNPs"
refFile="/1000G/GRCh38_EAS/1kGP_high_coverage_Illumina.EAS_unrelated.geno05.mind05.maf01.hwe1e6.autosome"
gctaAdd="/packages/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
${gctaAdd}  \
--bfile ${refFile} \
--cojo-file ${sumFile} \
--cojo-cond ${snpList} \
--out ${outFile} --threads 4