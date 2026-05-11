library(data.table)
rg <- fread("Rg_EA_53traits.txt")
index <- match(c("Type 2 diabetes","Refractive error","Schizophrenia","Prostate cancer","Glaucoma","Heart failure","Uterine fibroid","Body mass index","Triglycerides","Fasting glucose","Systolic blood pressure","HDL-C","Age at menarche","Drugs used in diabetes","Calcium channel blockers","Vitamin E intake"),rg$p2)
df_rg <- rg[index,]
df_rg$p2 <- factor(df_rg$p2,levels=df_rg$p2)

df_rg_EAS <- as.data.frame(df_rg)[,c(1,2,5,6:10)]
df_rg_EUR <- as.data.frame(df_rg)[,c(1,2,5,17:21)]
names(df_rg_EAS) <- c("p1","trait","Category","rg","se","z","p","padj")
df_rg_EAS$Ancestry <- "EAS"
df_rg_EAS$asterisk <- ""
df_rg_EAS$asterisk[df_rg_EAS$padj<0.05] <- "*"
names(df_rg_EUR) <- c("p1","trait","Category","rg","se","z","p","padj")
df_rg_EUR$Ancestry <- "EUR"
df_rg_EUR$asterisk <- ""
df_rg_EUR$asterisk[df_rg_EUR$padj<0.05] <- "*"
df_rg_all <- rbind(df_rg_EAS,df_rg_EUR)

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
library(ggplot2)
options(repr.plot.width = 7.2, repr.plot.height = 3.15)
df_rg_all$Ancestry <- factor(df_rg_all$Ancestry, levels = c("EAS", "EUR"))
p_rg <- ggplot(df_rg_all, aes(x = trait, y = rg, fill = Category, alpha = Ancestry)) +
  geom_col(width=0.7,position = position_dodge(width = 0.8), color = "white", linewidth = 0.1) +
  geom_errorbar(aes(ymin = rg - 1.96 * se, ymax = rg + 1.96 * se),width = 0.2, position = position_dodge(width = 0.8)) +
   geom_text(aes(
      y = ifelse(rg >= 0, 
                 rg + 1.96 * se + 0.02, 
                 rg - 1.96 * se - 0.02),
      label = asterisk,
      group = Ancestry
    ),
    position = position_dodge(width = 0.8),
    size = 8 / ggplot2::.pt,
    vjust = ifelse(df_rg_all$rg >= 0, 0, 1.5),
    alpha = 1,
    show.legend = FALSE
  )+
  scale_fill_manual(values = c(
  "Diseases"         = "#E64B35",
  "Biomarkers"       = "#E7B800",
  "Drug usage"       = "#4DBBD5",
  "Nutrition intake" = "#00A087"))+
  scale_alpha_manual(values = c("EAS" = 1.0, "EUR" = 0.4)) +
  ylab("Genetic correlation") +
  theme_nature_double()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_text(size = 8,face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white")
  ) + ylim(-0.45,0.7)

gsmr <- fread("GSMR_EA_53traits.gsmr")
gsmr$padj_EUR <- p.adjust(gsmr$p_EUR,method="fdr")
annot_sig <- function(x) {
if (is.na(x)) {
    sig <- ""
} else if (x > 0.05) {
    sig <- ""
} else if (x>0.01) {
    sig <- "*"
} else if (x> 0.001) {
    sig <- "**"
} else {
    sig <- "***"
};return(sig)}
gsmr_EAS <- as.data.frame(gsmr)[gsmr$Outcome %in% unique(gsub(" ","_",df_rg_all$trait)),c(1:7)]
names(gsmr_EAS) <- c("Outcome","Exposure","Beta","se","p","nsnp","padj")
gsmr_EAS$Exposure <- "EA (EAS)"
gsmr_EAS$asterisk <- sapply(gsmr_EAS$padj,annot_sig)
gsmr_EUR <- as.data.frame(gsmr)[gsmr$Outcome %in% unique(gsub(" ","_",df_rg_all$trait)),c(1,8:13)]
names(gsmr_EUR) <- c("Outcome","Exposure","Beta","se","p","nsnp","padj")
gsmr_EUR$Exposure <- "EA (EUR)"
gsmr_EUR$asterisk <- sapply(gsmr_EUR$padj,annot_sig)
df_gsmr <- rbind(gsmr_EAS,gsmr_EUR)
df_gsmr$trait <- gsub("_"," ",df_gsmr$Outcome)


library(ggplot2)
library(dplyr)
options(repr.plot.width = 7.2, repr.plot.height = 3.15)
df_gsmr$labels <-as.character(round(df_gsmr$Beta,digits=2))
df_gsmr$trait <- factor(df_gsmr$trait,levels = unique(df_rg_all$trait))
df_gsmr$Exposure <- factor(df_gsmr$Exposure,levels=c("EA (EUR)","EA (EAS)"))
p_gsmr2 <-  ggplot(df_gsmr,aes(y=Exposure,x=trait)) +
geom_tile(aes(fill=Beta),color = "white", lwd = 0.1, linetype = 1) +
geom_text(aes(label=asterisk),color = "black", size = 8 / ggplot2::.pt)+
scale_fill_gradient2(
    low = "#1583A1",
    high = "#C82512",
    mid = "#F2F2F2",      
    midpoint = 0, 
    na.value = "#888888",
    name = expression(beta[GSMR]))+
coord_equal()+
xlab("Outcome") +ylab("Exposure") +
theme_nature_double()+
theme(axis.text.x = element_text(size=7,angle=45,hjust=1),
    legend.position = "top",
    legend.key.width = unit(12, "mm"),
    legend.key.height = unit(4, "mm"),
    panel.border = element_blank())
p_gsmr2



library(patchwork)
options(repr.plot.width = 7.2, repr.plot.height = 6.3)
combined_plot <- p_rg/p_gsmr2 +
 plot_layout(nrow = 2, heights = c(2, 1))+
  plot_annotation(tag_levels = 'a') & 
  theme(
    plot.tag = element_text(size = 8, face = "bold", family = "Helvetica")
  )
combined_plot
















