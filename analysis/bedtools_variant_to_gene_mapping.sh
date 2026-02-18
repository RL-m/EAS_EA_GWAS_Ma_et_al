# bedtools_intersect.sh:
#module load bedtools/2.30.0
#file_a=$1
#file_b=$2
#file_c=$3
#bedtools \
#intersect -wa -wb \
#-a ${file_a} \
#-b ${file_b} > ${file_c}

snpFile=$1
outAdd=$2
prefix=$3
bedtools="bedtools_intersect.sh"

Gene_bed="gencode.v40.TSS.GRCh38.bed"
ABC_bed="ABC0.015_GRCh38_CHRALL.bed"
EpiMap_bed="EpiMap_links_by_group.ALL.CHRALL_GRCh38.bed"
Exon_bed="gencode.v40.GRCh38.exon.bed"
PCHiC_bed="PCHiC_Combined_GRCh38.bed"
RoadMap_bed="Roadmap_links_GRCh38.bed"


sh ${bedtools} ${snpFile} ${Exon_bed} ${outAdd}/${prefix}.V2G_Exon.txt

sh ${bedtools} ${snpFile} ${ABC_bed} ${outAdd}/${prefix}.V2G_ABC.txt

sh ${bedtools} ${snpFile} ${EpiMap_bed} ${outAdd}/${prefix}.V2G_EpiMap.txt

sh ${bedtools} ${snpFile} ${RoadMap_bed} ${outAdd}/${prefix}.V2G_RoadMap.txt

sh ${bedtools} ${snpFile} ${PCHiC_bed} ${outAdd}/${prefix}.V2G_PCHiC.txt