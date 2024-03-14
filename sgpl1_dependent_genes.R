#SGPL1 dependent genes
#Author: Dr Chris J Smith

#load packages
library(tidyverse)
library(matrixStats)
library(gplots)

#load data
genelist <- read.csv('skin_TPM_gene_names.csv')
ensembl_IDs <- read.csv('ENS_codes.csv')
control_vs_ko <- read.csv('skin_control_vs_ko_deseq2.csv')
control_vs_oe <- read.csv('skin_control_vs_oe_deseq2.csv')
tpm <- read.csv('skin_filtered_tpm_10.csv')

#makes genes rownames
genelist2 <- genelist[, -1:-2]
rownames(genelist2) <- genelist[, 1]

#Isolate genes that are of the right pattern
average_genotypes<- function(dataset) {
  control_mean <- apply(dataset[1:3], MARGIN = 1, mean)
  ko_mean <-  apply(dataset[4:6], MARGIN = 1, mean)
  oe_mean <- apply(dataset[7:9], MARGIN = 1, mean)
  means <- data.frame(control_mean, ko_mean, oe_mean)
  return(means)
}

average_of_genotypes <- average_genotypes(genelist2)

pattern_finder1 <- function(dataset) {
    ko_up_oe_down <- dataset %>%
    filter(ko_mean > control_mean & oe_mean < control_mean)
    ko_down_oe_up <- dataset %>%
    filter(ko_mean < control_mean & oe_mean > control_mean)
    return(ko_up_oe_down)
}

pattern_finder2 <- function(dataset) {
  ko_up_oe_down <- dataset %>%
    filter(ko_mean > control_mean & oe_mean < control_mean)
 ko_down_oe_up <- dataset %>%
    filter(ko_mean < control_mean & oe_mean > control_mean)
  return(ko_down_oe_up)
}

ko_up_oe_down <- pattern_finder1(average_of_genotypes)
ko_down_oe_up <- pattern_finder2(average_of_genotypes)

gene_pattern <- rbind(ko_up_oe_down, ko_down_oe_up)
gene_pattern <- cbind(rownames(gene_pattern), data.frame(gene_pattern, row.names=NULL))

names(gene_pattern)[1] <- "ENS"

gene_pattern <- inner_join(ensembl_IDs, gene_pattern, by = "ENS")

#Select genes that are logfold > 1
fold_change_genes <- control_vs_ko %>%
full_join(control_vs_oe, by = "ENS", suffix = c(".ko", ".oe")) %>%
mutate(abs_fc_ko = abs(log2FoldChange.ko)) %>%
mutate(abs_fc_oe = abs(log2FoldChange.oe)) %>%
filter(abs_fc_ko > 1  & abs_fc_oe > 1) %>%
select(ENS, Gene.ko) %>%
rename(Gene = Gene.ko)

#Select genes that are signifcant 
fold_change_genes <- control_vs_ko %>%
full_join(control_vs_oe, by = "ENS", suffix = c(".ko", ".oe")) %>%
filter(padj.ko < 0.05 | padj.oe < 0.05) %>%
select(ENS, Gene.ko) %>%
rename(Gene = Gene.ko)

#Find genes of the right pattern of the right fold change
sgpl1_dependent_genes <-  fold_change_genes %>%
semi_join(gene_pattern, by = "Gene")%>%
select(ENS, Gene)

tpm <- tpm %>%
rename(ENS = Geneid)

geneset_tpm <- semi_join(tpm, sgpl1_dependent_genes, by = "ENS")

#write csv
write.csv(geneset_tpm, 'SGPL1_dependent_genes_TPM_no_log_fold_ko_or_oe_skin.csv ')

#make heatmap
geneset_tpm_tibble<- geneset_tpm[,-1:-2]

rownames(geneset_tpm_tibble) <- geneset_tpm[,2]

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


heatmap_geneset_figure<-heatmap.2(z_score_matrix, dendrogram="row",Colv= "NA", scale="row", col=col, 
                          trace="none",ColSideColors=rep(c("blue","red","green"), each=3), density.info="density",
                          key = TRUE, keysize = 1.5, cexRow=0.5, cexCol=0.7,
                          main = "SGPL1 Dependent Genes" ,xlab = ,labRow = row.names(z_score_matrix))






