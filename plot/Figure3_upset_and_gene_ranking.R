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

library(data.table)
gene_score <- fread("EA_GWAS_EAS.gene_score.txt")
gene_df <- gene_score
names(gene_df)[13] <- "PoPS"
gene_df$MAGMA <- !is.na(gene_df$MAGMA)
gene_df$MBAT <- !is.na(gene_df$MBAT)
gene_df$V2G <- !is.na(gene_df$V2G)
gene_df$SMR <- !is.na(gene_df$SMR)
gene_df$COLOC <- !is.na(gene_df$COLOC)
gene_df$PoPS <- !is.na(gene_df$PoPS)

library(tidyverse)
comb_count <- gene_df %>%
mutate(
    COMB=pmap_chr(
    list(MAGMA,MBAT,V2G,SMR,COLOC,PoPS),
    \(lgl1,lgl2,lgl3,lgl4,lgl5,lgl6) {
        c("MAGMA","mBAT-combo","V2G","SMR","COLOC","PoPS")[c(lgl1,lgl2,lgl3,lgl4,lgl5,lgl6)]
    } %>% paste0(collapse = ",")
    )
) %>%
count(COMB) %>% 
mutate(
    COMB = fct_reorder(COMB,n,.desc = TRUE)
)
comb_count$score <- sapply(comb_count$COMB,function(x) {length(unlist(strsplit(as.character(x),",")))})
top_bar_colors <- c("#E0E0E0","#BDBDBD","#90CAF9","#42A5F5","#1E88E5","#0D47A1")
library(ggupset)
bar_top <- comb_count %>%
ggplot(aes(x=COMB,y=n,fill=score)) + geom_col() +
geom_text(aes(label = n),vjust = -0.2,color = "black",size = 5/.pt) +
scale_fill_gradientn(
    colors=top_bar_colors,
    name="No. of methods"
) +
scale_y_continuous(labels = NULL,expand=expansion(mult = c(0, 0.05)))+
theme_minimal(base_size=10) +
theme(
    panel.grid = element_blank(),
    axis.text.x=element_blank(),
    axis.title.y=element_text(color="black",size=8,face = "bold"),
    text = element_text(family = "Helvetica", size = 7),
    legend.position="top",
    legend.title = element_text(size = 6,face="bold"),
    legend.text = element_text(size = 5),
    legend.justification = c(1,-0.2),
    legend.key.height = unit(1, "mm"),
    legend.key.width = unit(4, "mm")
)+ labs(x="",y="Intersection size")

method_count <- data.frame(Method=c("V2G","MAGMA","mBAT-combo","SMR","COLOC","PoPS"),Count=apply(as.data.frame(comb_gene)[,c(7:10,12,13)],2,sum))
bar_left <- method_count %>%
mutate(Method=fct_reorder(Method,Count)) %>%
ggplot(aes(x=-Count,y=Method,fill=Method))+
geom_col()+
scale_fill_manual(values = c("COLOC"="#006064", "V2G"="#00838F","mBAT-combo"="#00ACC1", "PoPS"="#26C6DA", "SMR"="#66BB6A", "MAGMA"="#2E7D32"))+
geom_text(aes(label = Count),hjust=1,color="black",size=5/.pt)+
scale_x_continuous(labels=NULL,expand = expansion(mult=c(0.27,0))) +
theme_minimal(base_size=10) +
theme(
    panel.grid = element_blank(),
    axis.text.y=element_text(size=7,color="black",vjust=0.5,hjust=1),
    axis.text.x=element_blank(),
    text = element_text(family = "Helvetica", size = 7),
    margin(t = 5, r = 5, b = 5, l = 0, unit = "pt"),
    # Axis titles: Slightly larger (8pt) to make them "eury to tell"
    axis.title = element_text(size = 8, face = "bold", color = "black"),
    
    # Axis tick labels: Standard 7pt
    axis.text = element_text(size = 7, color = "black"),
    
    legend.position="none"
)+ labs(x="No. of genes",y="")


point_df <- comb_count %>% 
mutate(
    Method=map(
        COMB,
        ~str_split_1(as.character(.),",")
    )
) %>% 
unnest(Method) %>%
mutate(
    Method=factor(Method,levels=method_count %>% arrange(Count) %>% pull(Method)),
    COMB=factor(COMB,levels=levels(comb_count$COMB))
)
point_bottom <- point_df %>%
ggplot(aes(x=COMB,y=Method)) +
geom_line(aes(group=COMB),col="grey90") + 
geom_point(size=1,aes(color = Method))+
scale_color_manual(values = c("COLOC"="#006064", "V2G"="#00838F","mBAT-combo"="#00ACC1", "PoPS"="#26C6DA", "SMR"="#66BB6A", "MAGMA"="#2E7D32"))+
theme_minimal(base_size=10) +
theme(
    panel.grid = element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_blank(),
    axis.title.x=element_text(face="bold",size=8,color="black"),
    legend.position="none"
    
)+ labs(x="Intersection combinations",y="")

library(patchwork)
layout <- '
#AAAAAA
#AAAAAA
BCCCCCC'
f3 <- (bar_top) +
 (bar_left) +
  (point_bottom +
    scale_x_discrete(labels=NULL,expand=expansion(add=0.5)) +
        scale_y_discrete(expand=expansion(add=0.5))) +
   plot_layout(design=layout) &
   theme(
    panel.spacing = unit(1, "pt"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = -5, unit = "pt"))
options(repr.plot.width = 3.6, repr.plot.height = 3.15)
f3






