#Geneset to heatmap
#Author: Dr Chris J Smith

#load packages 
library(tidyverse)
library(matrixStats)
library(gplots)

#load data test
control_vs_ko <- read.csv('skin_control_vs_ko_deseq2.csv')
control_vs_oe <- read.csv('skin_control_vs_oe_deseq2.csv')
geneset <- read.csv('geneset.csv')
tpm <- read.csv('skin_TPM_gene_names.csv')

sig_ko_vs_control <- filter(control_vs_ko, padj <0.05)
sig_oe_vs_control <- filter(control_vs_oe, padj <0.05)

sig_geneset_ko <- inner_join(sig_ko_vs_control, geneset, by = "Gene")
sig_geneset_oe <- inner_join(sig_oe_vs_control, geneset, by = "Gene")

#write.csv(sig_geneset_ko, 'geneset_c_vs_KO.csv ')
#write.csv(sig_geneset_oe, 'geneset_c_vs_OE.csv ')

combined_geneset <- bind_rows(sig_geneset_ko, sig_geneset_oe)

geneset_tpm <- semi_join(tpm, combined_geneset, by = "ENS")

#write.csv(geneset_tpm, 'geneset._TPM.csv ')

#make variable 1 row names

geneset_tpm_tibble<- geneset_tpm[,-1:-2]

rownames(geneset_tpm_tibble) <- geneset_tpm[,2]

#import z score formula

create_z_score <- function(data.set) {
  row_mean <- rowMeans(as.matrix(data.set))
  row_sd <- rowSds(as.matrix(data.set))
  z_score <- (data.set - row_mean)/row_sd
  return(z_score)
}

#run z score formula

z_score <- create_z_score(geneset_tpm_tibble)

# set colour
col <- colorRampPalette(c("blue","white","red"))(30)
z_score_matrix <- as.matrix(z_score)

#plot heatmap
heatmap_geneset_figure<-heatmap.2(z_score_matrix, dendrogram="row",Colv= "NA", scale="row", col=col, 
                          trace="none",ColSideColors=rep(c("blue","red","green"), each=3), density.info="density",
                          key = TRUE, keysize = 1.5, cexRow=0.5, cexCol=0.7,
                          main = "KEGG Focal Adhesions" ,xlab = ,labRow = row.names(z_score_matrix))
