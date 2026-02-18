outAdd <- "/data"
dir.create(outAdd)
args=commandArgs(TRUE)
chr=as.numeric(args[1])
sumAdd=as.character(args[2])
qtlAdd=as.character(args[3])
qtlName=as.character(args[4])
probeAdd=as.character(args[5])
finalAdd=as.character(args[6])
qtlN=as.numeric(args[7])
index=as.numeric(args[8])
library(data.table)
library(coloc)
sumFile <- fread(sumAdd,header=T,col.names=c("SNP","A1","A2","AF1","BETA","SE","P","N"))
tmpName <- unlist(strsplit(sumAdd,"/"))
sumName <- tmpName[length(tmpName)]
### rm duplicated SNPs
dupSNP <- unique(sumFile$SNP[duplicated(sumFile$SNP)])
if (length(dupSNP)) {
sumFile <- sumFile[-which(sumFile$SNP %in% dupSNP),]
}
probeList <- fread(probeAdd,header=T)
#if (chr == 6) {
#mhcIndex <- which(probeList$Probe_bp <= 34000000 & probeList$Probe_bp >= 28000000)
#if (length(mhcIndex)) {
#probeList <- probeList[-mhcIndex,]
#}}
probeList <- probeList[which(probeList$Probe_Chr==chr),]
unit=ceiling(nrow(probeList)/10)
if (index<9) {
number=((unit*index+1):(unit*(index+1)))
} else {
number=((unit*index+1):(nrow(probeList)))
}
summaryList <- c()
dir.create(paste0(outAdd,"/query"))
dir.create(paste0(outAdd,"/result"))
for (probe in probeList$Probe[number]) {
system(paste0("
smr=/packages/smr_Linux
besdFile=",qtlAdd,"
probe=",probe,"
outFile=",outAdd,"/query/",qtlName,"_probe.",probe,"
${smr} --beqtl-summary ${besdFile} --probe ${probe} --query 1 --out ${outFile}"))
probeFile <- fread(paste0(outAdd,"/query/",qtlName,"_probe.",probe,".txt"),header=T)
gene <- unique(probeFile$Gene)
rmIndex <- which(probeFile$Freq == 1 | probeFile$Freq == 0)
if (length(rmIndex)) {
probeFile <- probeFile[-rmIndex,]
}
if (nrow(probeFile) <= 10) {
print(paste0("Not enough SNPs for the probe: ",qtlName,"_probe.",probe,"!"))
next
}
### keep common SNPs between probeFile & sumFile
commIndex <- match(probeFile$SNP,sumFile$SNP,nomatch=0)
if (length(which(commIndex != 0))<=10) {
print(paste0("Not enough SNPs for the probe: ",qtlName,"_probe.",probe,"!"))
next
}
commProbe <- probeFile[which(commIndex != 0),]
commSum <- sumFile[commIndex,]
### prepare input
list_sum <- list(beta=commSum$BETA,varbeta=(commSum$SE)^2,snp=commSum$SNP,position=commProbe$BP,type="quant",MAF=ifelse(commSum$AF1>1, 1-commSum$AF1, commSum$AF1),N=max(commSum$N))
# choose the corresponding command based on AF infomation in QTL file
list_xqtl <- list(beta=commProbe$b,varbeta=(commProbe$SE)^2,snp=commProbe$SNP,position=commProbe$BP,type="quant",MAF=ifelse(commProbe$Freq>1, 1-commProbe$Freq, commProbe$Freq),N=qtlN)
#list_xqtl <- list(beta=commProbe$b,varbeta=(commProbe$SE)^2,snp=commProbe$SNP,position=commProbe$BP,type="quant",MAF=list_sum[["MAF"]],N=qtlN)
check_dataset(list_sum)
check_dataset(list_xqtl)
### output results
coloc_res <- coloc.abf(dataset1=list_sum,dataset2=list_xqtl)
result <- data.frame(PROBE=probe,GENE=gene,CHR=probeFile$Probe_Chr[1],POS=probeFile$Probe_bp[1],PP0=coloc_res$summary[2],PP1=coloc_res$summary[3],PP2=coloc_res$summary[4],PP3=coloc_res$summary[5],PP4=coloc_res$summary[6],PP3_PP4=coloc_res$summary[5]+coloc_res$summary[6])
summaryList <- rbind(summaryList,result)
}
fwrite(summaryList,paste0(outAdd,"/result/",sumName,"_",qtlName,"_chr",chr,"_",index,".coloc.txt"),quote=F,sep="\t",col.names=T,row.names=F,na="NA")