# MED27 #
library(plotgardener)
library(plotgardenerData)
library(data.table)
outAdd <- ""
gwasFile <- fread("CARMA.EA_GWAS_EAS.txt.withoutAnnot.9.130978893_132978893.txt.gz",header=T)
prefix <- "9.130978893_132978893"
chr=9
START <- 131078893
END <- 132878893
eur <- fread("EA3.txt",header=T)
EAS_df <- data.frame(chrom=paste0("chr",gwasFile$CHR),pos=gwasFile$POS,snp=gwasFile$SNP,p=gwasFile$P,PIP=gwasFile$PIP,CS=gwasFile$CS)
fwrite(as.data.frame(EAS_df$snp),paste0(outAdd,"/LD.EAS.",prefix,".snplist"),quote=F,col.name=F,row.names=F,na="NA")
system(paste0("PLINK2=plink2
PLINK=/plink-1.9/plink
refFile=1kGP_high_coverage_Illumina.EAS_unrelated.geno05.mind05.maf01.hwe1e6.chr",chr,"
snpList=",outAdd,"/LD.EAS.",prefix,".snplist
ldFile=",outAdd,"/LD.EAS.",prefix,"
${PLINK2} --bfile ${refFile} --extract ${snpList} --make-bed --out ${ldFile}
${PLINK} --bfile ${ldFile} --r2 square spaces --keep-allele-order --threads 2 --out ${ldFile}
"))
EAS_bim <- fread(outAdd,"LD.EAS.",prefix,".bim",header=F,col.names=c("CHR","SNP","MOL","POS","A1","A2"))
EAS_ld <- as.matrix(fread(outAdd,"LD.EAS.",prefix,".ld",header=F))
colnames(EAS_ld) <- EAS_bim$SNP
rownames(EAS_ld) <- EAS_bim$SNP
EAS_top <- EAS_df[which.min(EAS_df$p),]
ld_r2 <- (EAS_ld[rownames(EAS_ld)==EAS_top$snp,])^2
index1 <- match(EAS_df$snp,names(ld_r2))
EAS_df$r2 <- ld_r2[index1]
EAS_df$ld <- "0.0-0.2"
EAS_df$col <- rgb(055,103,149,m=255)
EAS_df$ld[which(0.2 < EAS_df$r2&EAS_df$r2 <= 0.4)] <- "0.2-0.4"
EAS_df$col[which(0.2 < EAS_df$r2&EAS_df$r2 <= 0.4)] <- rgb(082,143,173,m=255)
EAS_df$ld[which(0.4 < EAS_df$r2&EAS_df$r2 <= 0.6)] <- "0.4-0.6"
EAS_df$col[which(0.4 < EAS_df$r2&EAS_df$r2 <= 0.6)] <- rgb(114,188,213,m=255)
EAS_df$ld[which(0.6 < EAS_df$r2&EAS_df$r2 <= 0.8)] <- "0.6-0.8"
EAS_df$col[which(0.6 < EAS_df$r2&EAS_df$r2 <= 0.8)] <- rgb(255,208,111,m=255)
EAS_df$ld[which(0.8 < EAS_df$r2&EAS_df$r2 <= 1)] <- "0.8-1.0"
EAS_df$col[which(0.8 < EAS_df$r2&EAS_df$r2 <= 1)] <- rgb(231,098,084,m=255)
EAS_df$y <- -log10(EAS_df$p)
EAS_df <- EAS_df[order(EAS_df$r2,decreasing=FALSE),]
# CARMA PIP #
pip_df <- EAS_df[,-4]
names(pip_df)[4] <- "p"
pip_df <- pip_df[pip_df$p > 0,]
pip_df <- pip_df[order(pip_df$p,decreasing=FALSE),]
# EUR GWAS #
eur <- eur[CHR==chr & POS <= END & POS >= START,]
EUR_df <- data.frame(chrom=paste0("chr",eur$CHR),pos=eur$POS,snp=eur$SNP,p=eur$P)
fwrite(as.data.frame(EUR_df$snp),paste0(outAdd,"/LD.EUR.",prefix,".snplist"),quote=F,col.name=F,row.names=F,na="NA")
system(paste0("PLINK2=plink2
PLINK=/plink-1.9/plink
refFile=BED_ukbEUR_imp_v3_INFO0.8_maf0.01_mind0.05_geno0.05_hwe1e6_10K_hg38_chr",chr,"
snpList=",outAdd,"/LD.EUR.",prefix,".snplist
ldFile=",outAdd,"/LD.EUR.",prefix,"
${PLINK2} --bfile ${refFile} --extract ${snpList} --make-bed --out ${ldFile}
${PLINK} --bfile ${ldFile} --r2 square spaces --keep-allele-order --threads 2 --out ${ldFile}
"))
EUR_bim <- fread(paste0(outAdd,"/LD.EUR.",prefix,".bim"),header=F,col.names=c("CHR","SNP","MOL","POS","A1","A2"))
EUR_ld <- as.matrix(fread(paste0(outAdd,"/LD.EUR.",prefix,".ld"),header=F))
colnames(EUR_ld) <- EUR_bim$SNP
rownames(EUR_ld) <- EUR_bim$SNP
EUR_top <- EUR_df[which.min(EUR_df$p),]
ld_r2 <- (EUR_ld[rownames(EUR_ld)==EAS_top$snp,])^2
index2 <- match(EUR_df$snp,names(ld_r2))
EUR_df$r2 <- ld_r2[index2]
EUR_df$ld <- "0.0-0.2"
EUR_df$col <- rgb(055,103,149,m=255)
EUR_df$ld[which(0.2 < EUR_df$r2&EUR_df$r2 <= 0.4)] <- "0.2-0.4"
EUR_df$col[which(0.2 < EUR_df$r2&EUR_df$r2 <= 0.4)] <- rgb(082,143,173,m=255)
EUR_df$ld[which(0.4 < EUR_df$r2&EUR_df$r2 <= 0.6)] <- "0.4-0.6"
EUR_df$col[which(0.4 < EUR_df$r2&EUR_df$r2 <= 0.6)] <- rgb(114,188,213,m=255)
EUR_df$ld[which(0.6 < EUR_df$r2&EUR_df$r2 <= 0.8)] <- "0.6-0.8"
EUR_df$col[which(0.6 < EUR_df$r2&EUR_df$r2 <= 0.8)] <- rgb(255,208,111,m=255)
EUR_df$ld[which(0.8 < EUR_df$r2&EUR_df$r2 <= 1)] <- "0.8-1.0"
EUR_df$col[which(0.8 < EUR_df$r2&EUR_df$r2 <= 1)] <- rgb(231,098,084,m=255)
EUR_df$y <- -log10(EUR_df$p)
EUR_df <- EUR_df[order(as.numeric(EUR_df$r2),decreasing=FALSE),]
# EAS eQTL #
EAS_x <- fread("JCTF_blood_eQTL.N1405.MED27.txt",header=T)
EAS_x <- EAS_x[BP>=START & BP<=END,]
EAS_x <- EAS_x[,c(2,3,1,14)]
names(EAS_x) <- c("chrom","pos","snp","p")
EAS_x$chrom <- paste0("chr",EAS_x$chrom)
index3 <- match(EAS_x$snp,EAS_df$snp,nomatch=0) 
EAS_x <- EAS_x[which(index3 != 0),]
EAS_x$col <- EAS_df$col[index3]
EAS_x$r2 <- EAS_df$r2[index3]
EAS_x <- EAS_x[order(EAS_x$r2,decreasing=FALSE),]
EAS_x$y <- -log10(EAS_x$p)
# annotation File #
abcFile <- fread("PCHiC_Combined_GRCh38.brain_ct.bed",header=F,col.names=c("CHR","START","END","GENE","Tissue","Feature"))
resultFile <- fread("EA_GWAS_EAS.V2G_PCHiC.txt",header=F)
annotFile <- abcFile[CHR==paste0("chr",chr) & GENE=="MED27",]
annotFile <- annotFile[annotFile$START >= START & annotFile$END <= END,]
annotFile$ID <- paste0(annotFile$CHR,":",annotFile$START,"_",annotFile$END)
annotFile <- annotFile[order(annotFile$Feature,decreasing=TRUE),]
annotFile <- annotFile[!duplicated(annotFile$ID),]
names(resultFile)[(ncol(resultFile)-5):ncol(resultFile)] <- names(abcFile)
resultFile <- resultFile[GENE=="MED27",]
highlightRegion <- resultFile[which.max(resultFile$Feature),]
geneList <- fread("gencode.v40.GRCh38.gene.annotation.bed",header=T)
geneList <- geneList[,c(1,6,8,7,3,4,5)]
names(geneList) <- c("CHR","PROBE","GENE","TYPE","GENESTART","GENEEND","ORIENTATION")
geneList <- geneList[which(geneList$TYPE=="protein_coding" & geneList$CHR %in% paste0("chr",1:22)),]
GENELIST <- geneList[which(geneList$CHR == paste0("chr",chr) & as.numeric(geneList$GENEEND) < END & as.numeric(geneList$GENESTART) > START),]
GENELIST$TSS <- ifelse(GENELIST$ORIENTATION=="+",GENELIST$GENESTART,GENELIST$GENEEND)
GENE <- GENELIST[GENE=="MED27",]
annot_df <- data.frame(chrom1=GENE$CHR,start1=GENE$TSS-500,end1=GENE$TSS+500,chrom2=annotFile$CHR,start2=annotFile$START,end2=annotFile$END)
annot_df$length <- abs(annot_df$start2-annot_df$start1)/1000
annot_df$h <- annot_df$length/max(annot_df$length)
highlight_df <- annot_df[which(annot_df$chrom2==highlightRegion$CHR & annot_df$start2==highlightRegion$START&annot_df$end2==highlightRegion$END),] 
index4 <- which(annot_df$start2==highlight_df$start2 & annot_df$end2==highlight_df$end2)
annot_df$col="grey"
annot_df$col[index4] <- "purple"
highlight_df <- annot_df[index4,]
annot_df <- rbind(as.data.frame(annot_df)[-index4,],highlight_df)
promoter <- data.frame(chrom="chr9",start=GENE$TSS-500,end=GENE$TSS+500)

pdf(paste0(outAdd,"/MED27.locus_plot.pdf"),width=3.6,height=6.3)
pageCreate(width = 9, height = 16, default.units = "cm")
pageGuideHide()
manhattanPlot1 <- plotManhattan(
  data = as.data.frame(EAS_df), chrom = paste0("chr",chr),
  chromstart = START,
  chromend = END,
  assembly = "hg38",
  fill = EAS_df$col,
  trans = "-log10",
  sigLine = TRUE, col = "grey",
  range=c(0,ceiling(max(EUR_df$y))+6),
  baseline.color="black",
  cex=0.2,
  lty = 2,
  leadSNP = list(
    snp = EAS_top$snp,
    pch=18,
    cex = 0.5,
    fill = "purple",
    fontsize = 0
  ),
  x = 1.5, y = 1.5, width = 6.5,
  height = 3,
  just = c("left", "top"),
  default.units = "cm"
)
annoYaxis(
  plot = manhattanPlot1,
  at = seq(0,ceiling(max(EUR_df$y))+2,by=round((max(EUR_df$y)+2-0)/4)),
  axisLine = TRUE, fontsize = 7
)
plotText(
    label = "EAS: EA", x = 1.6, y = 2, rot = 0,
    fontsize = 8, just = c("left","top"),fontface = "bold",
    default.units = "cm"
)
pipPlot1 <- plotManhattan(
  data = as.data.frame(pip_df), chrom = paste0("chr",chr),
  chromstart = START,
  chromend = END,
  assembly = "hg38",
  fill = pip_df$col,
  trans = "",
  sigLine = FALSE,
  range=c(0,max(pip_df$p)+0.07),
  baseline.color="black",
  cex=0.2,
  leadSNP = list(
    snp =EAS_top$snp,
    pch=18,
    cex = 0.5,
    fill = "purple",
    fontsize = 0
  ),
  x = 1.5, y = 4.5, width = 6.5,
  height = 1.5,
  just = c("left", "top"),
  default.units = "cm"
)
annoYaxis(
  plot = pipPlot1,
  at = c(0,round(max(pip_df$p),2)),
  axisLine = TRUE, fontsize = 7
)
manhattanPlot2 <- plotManhattan(
  data = as.data.frame(EUR_df), chrom = paste0("chr",chr),
  chromstart = START,
  chromend = END,
  assembly = "hg38",
  fill = EUR_df$col,
  trans = "-log10",
  sigLine = TRUE, col = "grey",
  range=c(0,ceiling(max(EUR_df$y))+6),
  baseline.color="black",
  cex=0.2,
  lty = 2,
  leadSNP = list(
    snp = EAS_top$snp,
    pch=18,
    cex = 0.5,
    fill = "purple",
    fontsize = 0
  ),
  x = 1.5, y = 6, width = 6.5,
  height = 3,
  just = c("left", "top"),
  default.units = "cm"
)
annoYaxis(
  plot = manhattanPlot2,
  at = seq(0,ceiling(max(EUR_df$y))+2,by=round((max(EUR_df$y)+2-0)/4)),
  axisLine = TRUE, fontsize = 7
)
plotText(
    label = "EUR: EA", x = 1.6, y = 6.5, rot = 0,
    fontsize = 8, just = c("left","top"),fontface = "bold",
    default.units = "cm"
)
manhattanPlot3 <- plotManhattan(
  data = as.data.frame(EAS_x), chrom = paste0("chr",chr),
  chromstart = START,
  chromend = END,
  assembly = "hg38",
  fill = EAS_x$col,
  trans = "-log10",
  sigLine = TRUE, col = "grey",
  range=c(0,ceiling(max(EAS_x$y))+20),
  baseline.color="black",
  cex=0.2,
  lty = 2,
  leadSNP = list(
    snp = EAS_top$snp,
    pch=18,
    cex = 0.5,
    fill = "purple",
    fontsize = 0
  ),
  x = 1.5, y = 9, width = 6.5,
  height = 3,
  just = c("left", "top"),
  default.units = "cm"
)
annoYaxis(
  plot = manhattanPlot3,
  at = seq(0,ceiling(max(EAS_x$y)+8)+1,by=round((max(EAS_x$y)+8)/4)),
  axisLine = TRUE, fontsize = 7
)
plotText(
    label = "EAS: eQTL\n(MED27)", x = 1.6, y = 9.5, rot = 0,
    fontsize = 8, just = c("left","top"),fontface = "bold",
    default.units = "cm"
)
annotPlot1 <- plotPairsArches(
    data = annot_df,
    chrom="chr9",
    chromstart = START,
    chromend = END,
    assembly = "hg38",
    fill = annot_df$col,linecolor = annot_df$col,
    flip=TRUE,clip=TRUE,
    archHeight = "h", alpha = 0.4,
    x = 1.5, y = 12, height = 1.5,width=6.5,
    just = c("left", "top"),
    default.units = "cm"
)
annotPlot3 <- plotRanges(
  data = promoter,
  chrom = "chr9", chromstart = START, chromend = END,
  assembly = "hg38",
  order = "random",collapse = TRUE,fill="red",
  x = 1.5, y = 12, width = 6.5, height = 0.5,
  just = c("left", "top"), default.units = "cm"
)
plotText(
    label = "PCHiC Interaction", x = 0.5, y = 12.5, rot = 0,
    fontsize = 8, just = c("left","top"),fontface = "bold",
    default.units = "cm"
)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
plotGenes(
  chrom = "chr9", chromstart = START, chromend = END,
  assembly = "hg38",
  geneHighlights = data.frame(
    "gene" = c("MED27"),
    "color" = c("#225EA8")
  ),
  x = 1.5, y =13.5, width = 6.5, height = 1,
  just = c("left", "top"), default.units = "cm",
  fontsize = 7
)
## Plot genome label
plotGenomeLabel(
  chrom = "chr9", chromstart = START, chromend = END,
  assembly = "hg38",
  x = 1.5, y = 14.5, length = 6.5, scale = "Mb",
  just = c("left", "top"), default.units = "cm",
  fontsize = 8
)
annoHighlight(
  plot = manhattanPlot1,
  chrom = "chr9",
  chromstart = EAS_top$pos-8000, chromend = EAS_top$pos+8000,
  y = 2, height = 12.5, just = c("top"),
  default.units = "cm",fill = "grey", alpha = 0.2,linecolor="grey"
)
plotText(
  label = EAS_top$snp, x = 3.25, y = 2.25,
  just = "left", default.units = "cm",fontsize = 7,
)
plotText(
    label = expression(-log[10]~'(p)'), x = 0.5, y = 3.25, rot = 90,
    fontsize = 8, just = "center",
    default.units = "cm"
)
plotText(
    label = "PIP", x = 0.5, y = 5.5, rot = 90,
    fontsize = 8, just = "center",
    default.units = "cm"
)
plotText(
    label = expression(-log[10]~'(p)'), x = 0.5, y = 7.75, rot = 90,
    fontsize = 8, just = "center",
    default.units = "cm"
)
plotText(
    label = expression(-log[10]~'(p)'), x = 0.5, y = 10.75, rot = 90,
    fontsize = 8, just = "center",
    default.units = "cm"
)
## Plot legend for LD scores
plotLegend(
  title=expression(LD ~ italic(r)^2),
  legend = c("Lead SNP","0.8-1.0","0.6-0.8","0.4-0.6","0.2-0.4","0-0.2"),
  fill = c("purple",rgb(231,098,084,m=255),rgb(255,208,111,m=255),rgb(114,188,213,m=255),rgb(082,143,173,m=255),rgb(055,103,149,m=255)), 
  cex = c(rep(0.5,5),1),
  pch = c(18, 19, 19, 19,19,19,19), border = FALSE, x = 7, y = 9.5,
  width = 1.5, height = 1.8, just = c("right", "top"),
  fontsize = 7,
  default.units = "cm"
)
dev.off()


# CBLN2 #
library(plotgardener)
library(plotgardenerData)
library(data.table)
outAdd <- ""
gwasFile <- fread("CARMA.EA_GWAS_EAS.txt.withoutAnnot.18.71543703_73543703.txt.gz",header=T)
prefix <- "18.71543703_73543703"
chr=18
START <- 71543703
END <- 73543703
eur <- fread("EA3.txt",header=T)
EAS_df <- data.frame(chrom=paste0("chr",gwasFile$CHR),pos=gwasFile$POS,snp=gwasFile$SNP,p=gwasFile$P,PIP=gwasFile$PIP,CS=gwasFile$CS)
fwrite(as.data.frame(EAS_df$snp),paste0(outAdd,"/LD.EAS.",prefix,".snplist"),quote=F,col.name=F,row.names=F,na="NA")
system(paste0("PLINK2=plink2
PLINK=/plink-1.9/plink
refFile=1kGP_high_coverage_Illumina.EAS_unrelated.geno05.mind05.maf01.hwe1e6.chr",chr,"
snpList=",outAdd,"/LD.EAS.",prefix,".snplist
ldFile=",outAdd,"/LD.EAS.",prefix,"
${PLINK2} --bfile ${refFile} --extract ${snpList} --make-bed --out ${ldFile}
${PLINK} --bfile ${ldFile} --r2 square spaces --keep-allele-order --threads 2 --out ${ldFile}
"))
EAS_bim <- fread(outAdd,"/LD.EAS.",prefix,".bim",header=F,col.names=c("CHR","SNP","MOL","POS","A1","A2"))
EAS_ld <- as.matrix(fread(outAdd,"/LD.EAS.",prefix,".ld",header=F))
colnames(EAS_ld) <- EAS_bim$SNP
rownames(EAS_ld) <- EAS_bim$SNP
EAS_top <- EAS_df[which.min(EAS_df$p),]
ld_r2 <- (EAS_ld[rownames(EAS_ld)==EAS_top$snp,])^2
index1 <- match(EAS_df$snp,names(ld_r2))
EAS_df$r2 <- ld_r2[index1]
EAS_df$ld <- "0.0-0.2"
EAS_df$col <- rgb(055,103,149,m=255)
EAS_df$ld[which(0.2 < EAS_df$r2&EAS_df$r2 <= 0.4)] <- "0.2-0.4"
EAS_df$col[which(0.2 < EAS_df$r2&EAS_df$r2 <= 0.4)] <- rgb(082,143,173,m=255)
EAS_df$ld[which(0.4 < EAS_df$r2&EAS_df$r2 <= 0.6)] <- "0.4-0.6"
EAS_df$col[which(0.4 < EAS_df$r2&EAS_df$r2 <= 0.6)] <- rgb(114,188,213,m=255)
EAS_df$ld[which(0.6 < EAS_df$r2&EAS_df$r2 <= 0.8)] <- "0.6-0.8"
EAS_df$col[which(0.6 < EAS_df$r2&EAS_df$r2 <= 0.8)] <- rgb(255,208,111,m=255)
EAS_df$ld[which(0.8 < EAS_df$r2&EAS_df$r2 <= 1)] <- "0.8-1.0"
EAS_df$col[which(0.8 < EAS_df$r2&EAS_df$r2 <= 1)] <- rgb(231,098,084,m=255)
EAS_df$y <- -log10(EAS_df$p)
EAS_df <- EAS_df[order(EAS_df$r2,decreasing=FALSE),]
# CARMA PIP #
pip_df <- EAS_df[,-4]
names(pip_df)[4] <- "p"
pip_df <- pip_df[pip_df$p > 0,]
pip_df <- pip_df[order(pip_df$p,decreasing=FALSE),]
# EUR GWAS #
eur <- eur[CHR==chr & POS <= END & POS >= START,]
EUR_df <- data.frame(chrom=paste0("chr",eur$CHR),pos=eur$POS,snp=eur$SNP,p=eur$P)
fwrite(as.data.frame(EUR_df$snp),paste0(outAdd,"/LD.EUR.",prefix,".snplist"),quote=F,col.name=F,row.names=F,na="NA")
system(paste0("PLINK2=plink2
PLINK=/plink-1.9/plink
refFile=BED_ukbEUR_imp_v3_INFO0.8_maf0.01_mind0.05_geno0.05_hwe1e6_10K_hg38_chr",chr,"
snpList=",outAdd,"/LD.EUR.",prefix,".snplist
ldFile=",outAdd,"/LD.EUR.",prefix,"
${PLINK2} --bfile ${refFile} --extract ${snpList} --make-bed --out ${ldFile}
${PLINK} --bfile ${ldFile} --r2 square spaces --keep-allele-order --threads 2 --out ${ldFile}
"))
EUR_bim <- fread(paste0(outAdd,"/LD.EUR.",prefix,".bim"),header=F,col.names=c("CHR","SNP","MOL","POS","A1","A2"))
EUR_ld <- as.matrix(fread(paste0(outAdd,"/LD.EUR.",prefix,".ld"),header=F))
colnames(EUR_ld) <- EUR_bim$SNP
rownames(EUR_ld) <- EUR_bim$SNP
EUR_top <- EUR_df[which.min(EUR_df$p),]
ld_r2 <- (EUR_ld[rownames(EUR_ld)==EAS_top$snp,])^2
index2 <- match(EUR_df$snp,names(ld_r2))
EUR_df$r2 <- ld_r2[index2]
EUR_df$ld <- "0.0-0.2"
EUR_df$col <- rgb(055,103,149,m=255)
EUR_df$ld[which(0.2 < EUR_df$r2&EUR_df$r2 <= 0.4)] <- "0.2-0.4"
EUR_df$col[which(0.2 < EUR_df$r2&EUR_df$r2 <= 0.4)] <- rgb(082,143,173,m=255)
EUR_df$ld[which(0.4 < EUR_df$r2&EUR_df$r2 <= 0.6)] <- "0.4-0.6"
EUR_df$col[which(0.4 < EUR_df$r2&EUR_df$r2 <= 0.6)] <- rgb(114,188,213,m=255)
EUR_df$ld[which(0.6 < EUR_df$r2&EUR_df$r2 <= 0.8)] <- "0.6-0.8"
EUR_df$col[which(0.6 < EUR_df$r2&EUR_df$r2 <= 0.8)] <- rgb(255,208,111,m=255)
EUR_df$ld[which(0.8 < EUR_df$r2&EUR_df$r2 <= 1)] <- "0.8-1.0"
EUR_df$col[which(0.8 < EUR_df$r2&EUR_df$r2 <= 1)] <- rgb(231,098,084,m=255)
EUR_df$y <- -log10(EUR_df$p)
EUR_df <- EUR_df[order(EUR_df$r2,decreasing=FALSE),]
# EAS mQTL #
EAS_x <- fread("EAS_mQTL_meta.CBLN2.cg25538235.txt",header=T)
EAS_x <- EAS_x[BP>=START & BP<=END,]
EAS_x <- EAS_x[,c(2,3,1,14)]
names(EAS_x) <- c("chrom","pos","snp","p")
EAS_x$chrom <- paste0("chr",EAS_x$chrom)
index3 <- match(EAS_x$snp,EAS_df$snp,nomatch=0) 
EAS_x <- EAS_x[which(index3 != 0),]
EAS_x$col <- EAS_df$col[index3]
EAS_x$r2 <- EAS_df$r2[index3]
EAS_x <- EAS_x[order(EAS_x$r2,decreasing=FALSE),]
EAS_x$y <- -log10(EAS_x$p)
# ABC annotation #
abcFile <- fread("ABC0.015_GRCh38_CHRALL.brain_ct.bed",header=F,col.names=c("CHR","START","END","GENE","Feature","Tissue"))
resultFile <- fread("EA_GWAS_EAS.V2G_ABC.txt",header=F)
annotFile <- abcFile[CHR==EAS_df$chrom[1] & GENE=="CBLN2",]
annotFile <- annotFile[annotFile$START > 71543703 & annotFile$END < 73543703,]
annotFile$ID <- paste0(annotFile$CHR,":",annotFile$START,"_",annotFile$END)
annotFile <- annotFile[order(annotFile$Feature,decreasing=TRUE),]
annotFile <- annotFile[!duplicated(annotFile$ID),]
names(resultFile)[(ncol(resultFile)-5):ncol(resultFile)] <- names(abcFile)
resultFile <- resultFile[GENE=="CBLN2",]
highlightRegion <- resultFile[which.max(resultFile$Feature),]
geneList <- fread("gencode.v40.GRCh38.gene.annotation.bed",header=T)
geneList <- geneList[,c(1,6,8,7,3,4,5)]
names(geneList) <- c("CHR","PROBE","GENE","TYPE","GENESTART","GENEEND","ORIENTATION")
geneList <- geneList[which(geneList$TYPE=="protein_coding" & geneList$CHR %in% paste0("chr",1:22)),]
GENELIST <- geneList[which(geneList$CHR == paste0("chr",chr) & as.numeric(geneList$GENEEND) < END & as.numeric(geneList$GENESTART) > START),]
GENELIST$TSS <- ifelse(GENELIST$ORIENTATION=="+",GENELIST$GENESTART,GENELIST$GENEEND)
GENE <- GENELIST[GENE=="CBLN2",]
annot_df <- data.frame(chrom1=GENE$CHR,start1=GENE$TSS-500,end1=GENE$TSS+500,chrom2=annotFile$CHR,start2=annotFile$START,end2=annotFile$END)
annot_df$length <- abs(annot_df$start2-annot_df$start1)/1000
annot_df$h <- annot_df$length/max(annot_df$length)
highlight_df <- annot_df[which(annot_df$chrom2==highlightRegion$CHR & annot_df$start2==highlightRegion$START&annot_df$end2==highlightRegion$END),] 
index4 <- which(annot_df$start2==highlight_df$start2 & annot_df$end2==highlight_df$end2)
annot_df$col="grey"
annot_df$col[index4] <- "purple"
highlight_df <- annot_df[index4,]
annot_df <- rbind(as.data.frame(annot_df)[-index4,],highlight_df)
promoter <- data.frame(chrom="chr18",start=GENE$TSS-500,end=GENE$TSS+500)

#svg(paste0(outAdd,"/CBLN2.locus_plot.svg"),width=3.6,height=6.3)
pdf(paste0(outAdd,"/CBLN2.locus_plot.pdf"),width=3.6,height=6.3)
#options(repr.plot.width = 3.6, repr.plot.height = 6.3)
pageCreate(width = 9, height = 16, default.units = "cm")
pageGuideHide()
#lotText(label = "A", fontsize = 26,
#            x = 1, y = 0, just = "left", default.units = "cm")
manhattanPlot1 <- plotManhattan(
  data = as.data.frame(EAS_df), chrom = paste0("chr",chr),
  chromstart = START,
  chromend = END,
  assembly = "hg38",
  fill = EAS_df$col,
  trans = "-log10",
  sigLine = TRUE, col = "grey",
  range=c(0,ceiling(max(EAS_df$y))+3.5),
  baseline.color="black",
  cex=0.2,
  lty = 2,
  leadSNP = list(
    snp = EAS_top$snp,
    pch=18,
    cex = 0.5,
    fill = "purple",
    fontsize = 0
  ),
  x = 1.5, y = 1.5, width = 6.5,
  height = 3,
  just = c("left", "top"),
  default.units = "cm"
)
annoYaxis(
  plot = manhattanPlot1,
  at = seq(0,ceiling(max(EAS_df$y))+1,by=round((max(EAS_df$y)-0)/4)),
  axisLine = TRUE, fontsize = 7
)
plotText(
    label = "EAS: EA", x = 1.6, y = 2, rot = 0,
    fontsize = 8, just = c("left","top"),fontface = "bold",
    default.units = "cm"
)
pipPlot1 <- plotManhattan(
  data = as.data.frame(pip_df), chrom = paste0("chr",chr),
  chromstart = START,
  chromend = END,
  assembly = "hg38",
  fill = pip_df$col,
  trans = "",
  sigLine = FALSE,
  range=c(0,max(pip_df$p)+0.07),
  baseline.color="black",
  cex=0.2,
  leadSNP = list(
    snp =EAS_top$snp,
    pch=18,
    cex = 0.5,
    fill = "purple",
    fontsize = 0
  ),
  x = 1.5, y = 4.5, width = 6.5,
  height = 1.5,
  just = c("left", "top"),
  default.units = "cm"
)
annoYaxis(
  plot = pipPlot1,
  at = c(0,round(max(pip_df$p),2)),
  axisLine = TRUE, fontsize = 7
)
manhattanPlot2 <- plotManhattan(
  data = as.data.frame(EUR_df), chrom = paste0("chr",chr),
  chromstart = START,
  chromend = END,
  assembly = "hg38",
  fill = EUR_df$col,
  trans = "-log10",
  sigLine = TRUE, col = "grey",
  range=c(0,ceiling(max(EAS_df$y))+3.5),
  baseline.color="black",
  cex=0.2,
  lty = 2,
  leadSNP = list(
    snp = EAS_top$snp,
    pch=18,
    cex = 0.5,
    fill = "purple",
    fontsize = 0
  ),
  x = 1.5, y = 6, width = 6.5,
  height = 3,
  just = c("left", "top"),
  default.units = "cm"
)
annoYaxis(
  plot = manhattanPlot2,
  at = seq(0,ceiling(max(EAS_df$y))+1,by=round((max(EAS_df$y)-0)/4)),
  axisLine = TRUE, fontsize = 7
)
plotText(
    label = "EUR: EA", x = 1.6, y = 6.5, rot = 0,
    fontsize = 8, just = c("left","top"),fontface = "bold",
    default.units = "cm"
)
manhattanPlot3 <- plotManhattan(
  data = as.data.frame(EAS_x), chrom = paste0("chr",chr),
  chromstart = START,
  chromend = END,
  assembly = "hg38",
  fill = EAS_x$col,
  trans = "-log10",
  sigLine = TRUE, col = "grey",
  range=c(0,ceiling(max(EAS_x$y))+24),
  baseline.color="black",
  cex=0.2,
  lty = 2,
  leadSNP = list(
    snp = EAS_top$snp,
    pch=18,
    cex = 0.5,
    fill = "purple",
    fontsize = 0
  ),
  x = 1.5, y = 9, width = 6.5,
  height = 3,
  just = c("left", "top"),
  default.units = "cm"
)
annoYaxis(
  plot = manhattanPlot3,
  at = seq(0,ceiling(max(EAS_x$y))+1,by=round((max(EAS_x$y)-0)/4)),
  axisLine = TRUE, fontsize = 7
)
plotText(
    label = "EAS: mQTL\n(cg25538235)", x = 1.6, y = 9.5, rot = 0,
    fontsize = 8, just = c("left","top"),fontface = "bold",
    default.units = "cm"
)
#params <- pgParams(
#    chrom="chr18",
#    chromstart=START,chromend=END,
#    assembly="hg38",
#    width=6.5
#)
annotPlot1 <- plotPairsArches(
    data = annot_df,
    chrom="chr18",
    chromstart = START,
    chromend = END,
    assembly = "hg38",
    fill = annot_df$col,linecolor = annot_df$col,
    flip=TRUE,clip=TRUE,
    archHeight = "h", alpha = 0.4,
    x = 1.5, y = 12, height = 1.5,width=6.5,
    just = c("left", "top"),
    default.units = "cm"
)
annotPlot3 <- plotRanges(
  data = promoter,
  chrom = "chr18", chromstart = START, chromend = END,
  assembly = "hg38",
  order = "random",collapse = TRUE,fill="red",
  x = 1.5, y = 12, width = 6.5, height = 0.5,
  just = c("left", "top"), default.units = "cm"
)
plotText(
    label = "ABC Interaction", x = 0.5, y = 12.5, rot = 0,
    fontsize = 8, just = c("left","top"),fontface = "bold",
    default.units = "cm"
)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
plotGenes(
  chrom = "chr18", chromstart = START, chromend = END,
  assembly = "hg38",
  geneHighlights = data.frame(
    "gene" = c("CBLN2"),
    "color" = c("#225EA8")
  ),
  x = 1.5, y =13.5, width = 6.5, height = 1,
  just = c("left", "top"), default.units = "cm",
  fontsize = 7
)
## Plot genome label
plotGenomeLabel(
  chrom = "chr18", chromstart = START, chromend = END,
  assembly = "hg38",
  x = 1.5, y = 14.5, length = 6.5, scale = "Mb",
  just = c("left", "top"), default.units = "cm",
  fontsize = 8
)
annoHighlight(
  plot = manhattanPlot1,
  chrom = "chr18",
  chromstart = EAS_top$pos-8000, chromend = EAS_top$pos+8000,
  y = 2, height = 12.5, just = c("top"),
  default.units = "cm",fill = "grey", alpha = 0.2,linecolor="grey"
)
plotText(
  label = EAS_top$snp, x = 3.25, y = 2.25,
  just = "left", default.units = "cm",fontsize = 7,
)
plotText(
    label = expression(-log[10]~'(p)'), x = 0.5, y = 3.25, rot = 90,
    fontsize = 8, just = "center",
    default.units = "cm"
)
plotText(
    label = "PIP", x = 0.5, y = 5.5, rot = 90,
    fontsize = 8, just = "center",
    default.units = "cm"
)
plotText(
    label = expression(-log[10]~'(p)'), x = 0.5, y = 7.75, rot = 90,
    fontsize = 8, just = "center",
    default.units = "cm"
)
plotText(
    label = expression(-log[10]~'(p)'), x = 0.5, y = 10.75, rot = 90,
    fontsize = 8, just = "center",
    default.units = "cm"
)
## Plot legend for LD scores
plotLegend(
  title=expression(LD ~ italic(r)^2),
  legend = c("Lead SNP","0.8-1.0","0.6-0.8","0.4-0.6","0.2-0.4","0-0.2"),
  fill = c("purple",rgb(231,098,084,m=255),rgb(255,208,111,m=255),rgb(114,188,213,m=255),rgb(082,143,173,m=255),rgb(055,103,149,m=255)), 
  cex = c(rep(0.5,5),1),
  pch = c(18, 19, 19, 19,19,19,19), border = FALSE, x = 7, y = 9.5,
  width = 1.5, height = 1.8, just = c("right", "top"),
  fontsize = 7,
  default.units = "cm"
)
dev.off()



