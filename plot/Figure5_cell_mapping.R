# plot ROSMAP scDRS result #
library(data.table)
library(ggplot2)
metaData <- as.data.frame(fread("ROSMAP_PFC.meta.csv",header=T))
scoreFile <- as.data.frame(fread("EA_GWAS_cross_ancestry.score.gz",header=T))
names(scoreFile)[1] <- "cell_name"
commIndex <- intersect(metaData$cell_name,scoreFile$cell_name)
scdrs_mama_df <- merge(metaData[match(commIndex,metaData$cell_name),],scoreFile[match(commIndex,scoreFile$cell_name),],by = "cell_name")

theme_nature_double <- function() {
  theme(
    text = element_text(family = "Helvetica", size = 7),
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    legend.title = element_text(size = 8,face="bold"),
    legend.text = element_text(size = 7),
    legend.background = element_blank(),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
   plot.margin = margin(0, 5, 5, 5, unit = "pt"),
   panel.grid = element_blank(),
   panel.background = element_rect(fill = "white"),
   plot.background = element_rect(fill = "white"),
  )
}

options(repr.plot.width = 3.6, repr.plot.height = 3.15)
scale_factor=1
p_umap <- ggplot(scdrs_mama_df) +
ggtitle('Cell type clustering\n(Human PFC)') + xlab("UMAP1") + ylab("UMAP2") +
geom_point(size=0.1,alpha=0.6,aes(x=umap1,y=umap2,col=ct)) + theme_nature_double()+
scale_color_manual(values = c("#999999", "#66C2A5","#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#E41A1C", "#377EB8","#4DAF4A", "#984EA3", "#FF7F00","#FFFF33","#A65628", "#F781BF")) +
theme(plot.title = element_text(size = 8, face = "bold", family = "Helvetica",hjust = 0.5, vjust = 0.5),
            legend.position.inside=c(0,0),
            legend.direction = "vertical",
            legend.justification = c(0,0),
            legend.key.width = unit(2, "mm"),
            legend.key.height = unit(3, "mm"),
            ) + labs(color = "Cell type") +guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

options(repr.plot.width = 3.6, repr.plot.height = 3.15)
color_pval = rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695"))
scdrs_mama_df <- scdrs_mama_df[order(scdrs_mama_df$nlog10_pval,decreasing=FALSE),]
p_assoc <- ggplot(scdrs_mama_df) +
ggtitle('EA association\n(Multi-ancestry meta-analysis)')+
geom_point(size=0.1,alpha=0.6,aes(x=umap1,y=umap2,col=nlog10_pval))+
scale_color_gradientn(colours = color_pval,
                            na.value = '#636363',
                            labels = scales::number_format(accuracy = 1)) +
guides(color = guide_colorbar(title.position = "top")) + theme_nature_double()+
theme(plot.title = element_text(size = 8, face = "bold", family = "Helvetica",hjust = 0.5, vjust = 0.5),
            legend.position.inside=c(0,0),
            legend.direction = "vertical",
            legend.justification = c(0, 0),
            legend.key.width = unit(2, "mm"),
            legend.key.height = unit(5, "mm"),
            )+labs(x="UMAP1",y="UMAP2",color = expression(-log[10]~'(p)'))


scdrs_ct <- fread("EA_GWAS_cross_ancestry.scDRS.summary.txt")
scdrs_ct$padj <- p.adjust(scdrs_ct$assoc_mcp,method="fdr")
scdrs_ct$logp <- -log10(scdrs_ct$padj)
scdrs_ct <- scdrs_ct[order(scdrs_ct$assoc_mcz,decreasing=FALSE),]
scdrs_ct$group <- factor(scdrs_ct$group,levels=scdrs_ct$group)
options(repr.plot.width = 3.2, repr.plot.height = 3.15)
fig_scdrs_ct <- ggplot(scdrs_ct) + xlab(expression(-log[10]~'(p)')) + ylab("Cell type (Human PFC)")+
    geom_bar( aes(x=logp, y=group), stat="identity", fill="#D81B60", alpha=0.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey")+ 
    theme_nature_double()


library(data.table)
region_ct_assoc <- fread("EA_GWAS_cross_ancestry.rg_ct_cauchy.txt")
region_ct_assoc <- region_ct_assoc[order(region_ct_assoc$p,decreasing=FALSE),]
region_ct_assoc <- region_ct_assoc[region_ct_assoc$CT != "OPC",]
region_ct_assoc$cell_type <- factor(region_ct_assoc$cell_type,levels=unique(region_ct_assoc$cell_type))
region_ct_assoc$region <- factor(region_ct_assoc$region,levels=rev(unique(region_ct_assoc$region)))
region_ct_assoc$CT <- factor(region_ct_assoc$CT,levels=unique(region_ct_assoc$CT))
region_ct_assoc$rg <- factor(region_ct_assoc$rg,levels=rev(unique(region_ct_assoc$rg)))
library(ggplot2)
options(repr.plot.width = 5, repr.plot.height = 3.15)
library(stringr)
color_pval = c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")
region_ct_heatmap <- ggplot(region_ct_assoc, aes(x = CT, y = rg, fill = logp)) +
  geom_tile(color = "white", lwd = 0.5, linetype = 1) +
  scale_y_discrete(position="right")+
  geom_text(size=7,aes(label = asterisk), 
            color = "white",
            fontface = "bold",vjust=0.8)+
  scale_fill_gradientn(colors=rev(color_pval),name = expression(-log[10]~'(p)'),na.value = "#bababa") +
  theme_nature_double()+
  theme(plot.margin =margin(c(t=5,r=1,b=5,l=10), unit="pt"),
    axis.text.x = element_text(angle = 30,hjust = 1),
  axis.text.y = element_text(angle = 30,vjust =0.1),
  legend.position = 'right',
  legend.direction = "vertical",
  legend.justification = c(0,0),
  legend.key.width = unit(2, "mm"),
  legend.key.height = unit(5, "mm")) + labs(x="Cell type (Mouse brain)",y="Brain region (Mouse brain)")

# plot Stereo-seq gsMap result #
library(data.table)
mamaResult <- "Adult_Mouse_Brain_EA_GWAS_cross_ancestry.csv.gz"
obs_meta_MAMA = readRDS('Adult_Mouse_brain_cell_bin_spatail_ldsc.rds')
mamaEA = as.data.frame(fread(mamaResult))
rownames(mamaEA) = mamaEA$spot
mamaEA$p.adj <- p.adjust(mamaEA$p,method="fdr")
rownames(obs_meta_MAMA) <- obs_meta_MAMA$cell_name
obs_meta_MAMA = obs_meta_MAMA[,c('cell_name','annotation','region','x','y')]
obs_meta_MAMA$cell_type=obs_meta_MAMA$annotation
obs_meta_MAMA$EA_MAMA <- -log10(mamaEA$p[match(rownames(obs_meta_MAMA),rownames(mamaEA))])
obs_meta_MAMA$size = 'A'
obs_meta_MAMA$size[obs_meta_MAMA$EA_MAMA>6 & obs_meta_MAMA$EA_MAMA<8] = 'B'
obs_meta_MAMA$size[obs_meta_MAMA$EA_MAMA>8] = 'C'
obs_meta_MAMA <- obs_meta_MAMA[complete.cases(obs_meta_MAMA),]
obs_meta_MAMA <- obs_meta_MAMA[obs_meta_MAMA$annotation!="Unknown",]
library(ggplot2)
options(repr.plot.width = 3.6, repr.plot.height = 3.15)
scale_factor=1
color_pval = c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")
color_pval = rev(color_pval)
fig_mouse_brain_ST = ggplot(obs_meta_MAMA,aes(x=x,y=y))+
    ggtitle('Brain region\n(Mouse brain)')+labs(color = "Brain region")+
    geom_point(aes(col=region,size=size))+
    scale_size_manual(values = c(0.01,0.05,0.1),guide = 'none')+
    scale_color_manual(values = c("#999999", "#66C2A5","#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#E41A1C", "#377EB8","#4DAF4A", "#984EA3", "#FF7F00","#FFFF33","#A65628", "#F781BF"))+
    theme_nature_double()+
    theme(plot.title = element_text(size = 8, face = "bold", family = "Helvetica",hjust = 0.5, vjust = 0.5),
            axis.line=element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.5,0.8,0.2,0), "cm"),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = 'right',
            legend.direction = "vertical",
            legend.justification = c(0, 0),
            legend.key.width = unit(2, "mm"),
            legend.key.height = unit(3, "mm"),
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"))+guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 1)))

obs_meta_MAMA <- obs_meta_MAMA[order(obs_meta_MAMA$EA_MAMA,decreasing=FALSE),]
fig_EA_MAMA_gsmap = ggplot(obs_meta_MAMA,aes(x=x,y=y))+
    ggtitle('EA association\n(Multi-ancestry meta-analysis)')+
    geom_point(aes(col=EA_MAMA,size=size))+
    scale_size_manual(values = c(0.01,0.05,0.1),guide = 'none')+
    scale_color_gradientn(colours = color_pval,
                            na.value = '#636363',
                            labels = scales::number_format(accuracy = 1))+
    guides(color = guide_colorbar(title.position = "top"))+
    theme_nature_double()+
    theme(plot.title = element_text(size = 8, face = "bold", family = "Helvetica",hjust = 0.5, vjust = 0.5),
            axis.line=element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.5,0.8,0.2,0), "cm"),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = 'right',
            legend.direction = "vertical",
            legend.justification = c(0,0),
            legend.key.width = unit(2, "mm"),
            legend.key.height = unit(5, "mm"),
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"))+
            labs(color = expression(-log[10]~'(p)'))

p1 <- p_umap
p2 <- p_assoc
p3 <- fig_scdrs_ct
p4 <- region_ct_heatmap
p5 <- fig_mouse_brain_ST
p6 <- fig_EA_MAMA_gsmap
library(patchwork)
options(repr.plot.width = 7.2, repr.plot.height = 9)
row1 <- p1 + p2
row2 <- p3 + p4 + plot_layout(widths = c(1, 1.2))
row3 <- p5 + p6
final_plot <- (row1 / row2 / row3) + 
  plot_annotation(tag_levels = 'a') & 
  theme(    plot.margin = margin(t=2, r=0, b=5, l=7,unit="pt"),
    plot.tag = element_text(size = 8, face = "bold", family = "Helvetica")
  )
final_plot








