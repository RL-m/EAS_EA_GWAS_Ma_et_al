# Figure 2 Manhattan plot
library(ggplot2)
theme_nature_double <- function() {
  #theme_classic() + # Removes grey background and grids
  theme(
    # Use Helvetica or Arial as preferred by high-impact journals
    text = element_text(family = "Helvetica", size = 7),
    
    # Axis titles: Slightly larger (8pt) to make them "easy to tell"
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    
    # Axis tick labels: Standard 7pt
    axis.text = element_text(size = 7, color = "black"),
    
    # Legend: Title at 8pt, keys at 7pt
    legend.title = element_text(size = 8,face="bold"),
    legend.text = element_text(size = 7),
    legend.background = element_blank(),
    legend.key.size = unit(4, "mm"),
    
    # Lines: Clean axis lines
    axis.line = element_line(linewidth = 0.3), # standard ~0.85pt
    axis.ticks = element_line(linewidth = 0.3),
   plot.margin = margin(0, 5, 5, 5, unit = "pt"),    # margin(top, right, bottom, left)
   panel.grid = element_blank()
  )
}
# p1: EAS GWAS with EAS-specific associations #
library(data.table)
sumFile <- fread("EA_GWAS_EAS.txt",header=T)
plotFile <- sumFile[,c(1,2,3,7,9)]
clump <- fread("EA_GWAS_EAS.after_cond_EA4_sigSNPs.cma.cojo.clumped")
highlight_all <- c()
for (i in 1:nrow(clump)) {
indexSNP <- sumFile[which(sumFile$SNP==clump$SNP[i]),]
highlight <- sumFile[sumFile$CHR==indexSNP$CHR & sumFile$POS<=indexSNP$POS+1000000 & sumFile$POS>=indexSNP$POS-1000000,]
highlight_all <- rbind(highlight_all,highlight)
}
sumFile$CHR <- paste0("chr",sumFile$CHR)
highlight_all$CHR <- paste0("chr",highlight_all$CHR)

library(topr)
library(ggplot2)
library(scales)
options(repr.plot.width = 7.2, repr.plot.height = 3.15)
p1_pre <- manhattan(list(sumFile, highlight_all),
          color = c("#003f5c","#ffa600"), 
          even_no_chr_lightness = c(0.7,0.5), 
          size=c(0.2,0.3),
          sign_thresh = 5e-08,
          sign_thresh_color = "black",
          sign_thresh_size=0.3,
          sign_thresh_label_size=2.5,
          label_size = 3,
          verbose = FALSE,
          show_legend=FALSE
          ) + 
  theme_minimal() +
  theme(panel.grid = element_blank(),
    axis.line = element_line(color = "grey"),
        legend.position = "top") +
  theme_nature_double() + scale_y_continuous(
  # Force breaks to be at integers
  breaks = breaks_pretty(), 
  # Round labels to the nearest whole number
  labels = function(x) ceiling(x) 
)
p1_pro <- p1_pre+ geom_point(aes(x = 0, y = 0, fill = "EAS-specific associations"), 
             shape = 21, size = 0, color = NA) + 
  scale_fill_manual(name = NULL, values = c("EAS-specific associations" = "#ffa600")) +
  guides(fill = guide_legend(override.aes = list(size = 1.5, shape = 21))) + 
  theme(
    legend.position = "top", 
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 7, family = "Helvetica"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.margin = margin(t = 0, b = 0, l = 0),
    plot.margin = margin(t = 3, r = 10, b = 5, l = 5)
  )
p1_pro


# p2: EA GWAS cross-ancestry with novel loci #
library(data.table)
metaFile <- fread("EA_GWAS_cross_ancestry.txt")
metaClump <- fread("EA_GWAS_cross_ancestry.clumped")
overlap <- fread("EA_GWAS_cross_ancestry_EA4_intersect.txt",header=F)
newClump <- metaClump[!metaClump$SNP %in% unique(overlap$V4),]

highlight_meta <- c()
for (i in 1:nrow(newClump)) {
indexSNP <- metaFile[which(metaFile$SNP==newClump$SNP[i]),]
highlight <- metaFile[metaFile$CHR==indexSNP$CHR & metaFile$POS<=indexSNP$POS+1000000 & metaFile$POS>=indexSNP$POS-1000000,]
highlight_meta <- rbind(highlight_meta,highlight)
}
metaFile$CHR <- paste0("chr",metaFile$CHR)
highlight_meta$CHR <- paste0("chr",highlight_meta$CHR)

library(topr)
library(ggplot2)
library(scales)
options(repr.plot.width = 7.2, repr.plot.height = 3.15)
p2_pre <- manhattan(list(metaFile, highlight_meta),
          color = c("#003f5c","#ffa600"), 
          even_no_chr_lightness = c(0.7,0.5), 
          size=c(0.2,0.3),
          sign_thresh = 5e-08,
          sign_thresh_color = "black",
          sign_thresh_size=0.3,
          sign_thresh_label_size=2.5,
          label_size = 3,
          verbose = FALSE,
          show_legend=FALSE
          ) + 
  theme_minimal() +
  theme(panel.grid = element_blank(),
    axis.line = element_line(color = "grey"),
        legend.position = "top") +
  theme_nature_double() + scale_y_continuous(
  breaks = breaks_pretty(), 
  labels = function(x) ceiling(x) 
)
p2_pro <- p2_pre + geom_point(aes(x = 0, y = 0, fill = "Novel loci compared with Okbay et. al."), 
             shape = 21, size = 0, color = NA) + 
  scale_fill_manual(name = NULL, values = c("Novel loci compared with Okbay et. al." = "#ffa600")) +
  guides(fill = guide_legend(override.aes = list(size = 1.5, shape = 21))) +
  theme(
    legend.position = "top", 
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 7, family = "Helvetica"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.margin = margin(t = 0, b = 0, l = 0),
    plot.margin = margin(t = 3, r = 10, b = 5, l = 5)
  )
p2_pro


library(ggpubr)
options(repr.plot.width = 7.2, repr.plot.height = 6.3)
f2_pro <- ggarrange(p1_pro,p2_pro,nrow=2,ncol=1,heights=c(1,1),labels=c("a","b"),font.label = list(size = 8, color = "black", face = "bold"))
f2_pro










