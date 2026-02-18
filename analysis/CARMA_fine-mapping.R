args=commandArgs(TRUE)
chr=as.numeric(args[1])
POS=as.numeric(args[2])
mainAdd=args[3]
sumFileName=args[4]
trait_pre <- unlist(strsplit(sumFileName,"/"))
trait <- trait_pre[length(trait_pre)]
START <- POS-1000000
END <- POS+1000000
library(data.table)
rawRef <- fread(paste0("/LD_reference/1000G/1KGP3_3202_GRCh38/EAS/1kGP_high_coverage_Illumina.EAS_unrelated.geno05.mind05.maf01.hwe1e6.chr",chr,".bim"),header=F,col.names=c("CHR","SNP","MOL","POS","A1","A2"))
ldAdd <- "/data/LD"
dir.create(ldAdd)
rawSum <- fread(sumFileName,header=T)
rawSum <- rawSum[which(rawSum$CHR == chr),]
rawSum$Z <- rawSum$BETA/rawSum$SE
## format the inputs
prefix <- paste0(trait,".withoutAnnot.",chr,".",START,"_",END)
sumFile <- rawSum[which(rawSum$POS >= START & rawSum$POS <= END),]
refFile <- rawRef[which(rawRef$POS >= START & rawRef$POS <= END),]
commSNP <- intersect(sumFile$SNP,refFile$SNP)
if (length(commSNP) <= 10) {
print(paste0("Not enough SNP in region: ",prefix,"!"))
next
}
### Keep the SNP order as refFile
commSum <- sumFile[match(commSNP,sumFile$SNP),]
commRef <- refFile[match(commSNP,refFile$SNP),]
## Flip alleles to make consistent sumFile & LD
flipIndex <- which(commSum$A1 != commRef$A1)
A0 <- commSum$A1[flipIndex]
commSum$A1[flipIndex] <- commSum$A2[flipIndex]
commSum$A2[flipIndex] <- A0
commSum$Z[flipIndex] <- -commSum$Z[flipIndex]
commSum$BETA[flipIndex] <- -commSum$BETA[flipIndex]
commSum$AF1[flipIndex] <- 1-commSum$AF1[flipIndex]
keepIndex <- which(commSum$A1 == commRef$A1 & commSum$A2 == commRef$A2)
commSum <- commSum[keepIndex,]
commRef <- commRef[keepIndex,]
fwrite(as.data.frame(commRef$SNP),paste0(ldAdd,"/LD.",prefix,".snplist"),quote=F,col.names=F,row.names=F,na="NA")
system(paste0("Plink19Add=/packages/plink-1.9/plink
refFile=/LD_reference/1000G/1KGP3_3202_GRCh38/EAS/1kGP_high_coverage_Illumina.EAS_unrelated.geno05.mind05.maf01.hwe1e6.chr",chr,"
snpList=",ldAdd,"/LD.",prefix,".snplist
ldFile=",ldAdd,"/LD.",prefix,"
${Plink19Add} --bfile ${refFile} --extract ${snpList} --keep-allele-order --make-bed --out ${ldFile}
${Plink19Add} --bfile ${ldFile} --r square spaces --keep-allele-order --threads 2 --out ${ldFile}"))
## Get final inputs
bimFile <- fread(paste0(ldAdd,"/LD.",prefix,".bim"),header=F,col.names=c("CHR","SNP","MOL","POS","A1","A2"))
ld <- as.matrix(fread(paste0(ldAdd,"/LD.",prefix,".ld"),header=F))
finalIndex <- which(complete.cases(ld))
if (length(finalIndex) <= 10) {
print(paste0("Not enough SNP in region: ",prefix,"!"))
next
}
ld <- ld[finalIndex,finalIndex]
refFinal <- bimFile[finalIndex,]
sumFinal <- commSum[match(refFinal$SNP,commSum$SNP),]
## Run CARMA
library(magrittr)
library(dplyr)
library(devtools)
library(R.utils)
library(CARMA)
library(data.table)
z.list<-list()
ld.list<-list()
lambda.list<-list()
annot.list<-list()
z.list[[1]]<-sumFinal$Z
ld.list[[1]]<-ld
lambda.list[[1]]<-1
CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,outlier.switch=T)
sumstat.result = sumFinal %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
}}
warnings()
fwrite(sumstat.result,paste0(mainAdd,"/CARMA.",prefix,".txt.gz"), sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")
file.remove("/data/LD")