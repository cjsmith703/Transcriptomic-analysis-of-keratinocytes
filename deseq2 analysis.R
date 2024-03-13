#DESEQ2 for Bulk RNAseq of SGPL1 KO keratinocytes
#Author: Dr Chris J Smith

#load packages
library(DESeq2)
library(GEOquery)
library(canvasXpress)
library(tidyverse)
library(clinfun)
library(GGally)
library(factoextra)
library(ashr)

#download clinical data from GEO
data <- getGEO(GEO = "GSE207499")
head(data)

clindata <- data[["GSE207499_series_matrix.txt.gz"]]@phenoData@data

#Alternatively read csv of clinical data if not uploaded to GEO repository

#read raw counts file
raw_counts <- read.csv("skin_feature_counts.csv")
#add ensembl codes
ens <- read.csv("ENS_codes.csv")


#Rename rows & columns so that row names in clindata = colnames in raw_counts
rownames(clindata) <- clindata$title

colnames(raw_counts)[c(2:10)] <- clindata$title

rownames(raw_counts) <- raw_counts$Geneid
raw_counts <- raw_counts[, -1]

raw_counts <- as.matrix(raw_counts)

#Check rows = columns, both should be TRUE
all(rownames(clindata) %in% colnames(raw_counts))

all(colnames(raw_counts) %in% rownames(clindata))

#set variables as factor and rename them
clindata$genotype <- clindata$`genotype:ch1`

clindata$genotype[clindata$genotype == "WT"] <- "SGPL1_Control"
clindata$genotype[clindata$genotype == "SGPL1 KO"] <- "SGPL1_KO"
clindata$genotype[clindata$genotype == "SGPl1 OE"] <- "SGPL1_OE"

clindata$genotype <- as.factor(clindata$genotype)

class(clindata$genotype)

#construct deseq2 data set (dds)
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = clindata,
                              design = ~genotype)

head(dds)

#filter for counts > 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Perform deseq2
dds <- DESeq(dds)

levels(dds$genotype)

res <- results(dds, contrast = c("genotype", "SGPL1_Control", "SGPL1_KO"))
res

#Order by p value
resOrdered <- res[order(res$pvalue), ]
resOrdered

#How many significant genes
sum(res$padj < 0.05, na.rm = TRUE)

resultsNames(dds)

#Log fold shrinkage - common consensus to do this way
control_vs_ko <- lfcShrink(dds,
                           contrast = c("genotype", "SGPL1_KO", "SGPL1_Control"),
                           type = "ashr")
control_vs_ko

control_vs_oe <- lfcShrink(dds,
                           contrast = c("genotype", "SGPL1_OE", "SGPL1_Control"),
                           type = "ashr")
control_vs_oe

ko_vs_oe <- lfcShrink(dds,
                      contrast = c("genotype", "SGPL1_OE", "SGPL1_KO"),
                      type = "ashr")
ko_vs_oe

#Convert to data fame and add gene names
control_vs_ko_deseq2 <- control_vs_ko %>%
  as.data.frame() %>%
  arrange(padj) %>%
  tibble::rownames_to_column() %>%
  rename(ENS = rowname) %>%
  inner_join(ens, by = "ENS") %>%
  relocate(Gene, .after = "ENS")

#write.csv(control_vs_ko_deseq2, 'skin_control_vs_ko_deseq2.csv ')
  
control_vs_oe_deseq2 <- control_vs_oe %>%
  as.data.frame() %>%
  arrange(padj) %>%
  tibble::rownames_to_column() %>%
  rename(ENS = rowname) %>%
  inner_join(ens, by = "ENS") %>%
  relocate(Gene, .after = "ENS")

#write.csv(control_vs_oe_deseq2, 'skin_control_vs_oe_deseq2.csv ')

ko_vs_oe_deseq2 <- ko_vs_oe %>%
  as.data.frame() %>%
  arrange(padj) %>%
  tibble::rownames_to_column() %>%
  rename(ENS=rowname) %>%
  inner_join(ens, by = "ENS") %>%
  relocate(Gene, .after = "ENS")

#write.csv(ko_vs_oe_deseq2, 'skin_ko_vs_oe_deseq2.csv ')


#How many significant genes
sum(control_vs_ko$padj < 0.05, na.rm = TRUE)
sum(control_vs_oe$padj < 0.05, na.rm = TRUE)
sum(ko_vs_oe$padj < 0.05, na.rm = TRUE)

#Info on columns
mcols(control_vs_ko)$description

##Plots

#MA-plots

plotMA(control_vs_ko, ylim = c(-6, 6))

plotMA(control_vs_oe, ylim = c(-6, 6))

plotMA(ko_vs_oe, ylim = c(-6, 6))


#Identifies rows by clicking on the plot
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]


#PCA Plot
colours <- c(SGPL1_Control = "blue",
             SGPL1_KO = "red",
             SGPL1_OE = "green")
vsdata <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsdata, intgroup ="genotype", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
  scale_colour_manual(values = colours)+
  geom_point(size = 2) +
  xlab(paste0("PC1: ", percentVar[1],"% variance")) +
  ylab(paste0("PC2: ", percentVar[2],"% variance")) +
  coord_fixed() +
  theme_minimal()

#Plot counts - plots the counts for gene with the lowest p value
plotCounts(dds, gene = which.min(control_vs_ko$padj), intgroup = "genotype")

plotCounts(dds, gene = which.min(control_vs_oe$padj), intgroup = "genotype")

plotCounts(dds, gene = which.min(ko_vs_oe$padj), intgroup ="genotype")


cko_plot <- plotCounts(dds, gene = which.min(control_vs_ko$padj),
                      intgroup = "genotype", returnData = TRUE)
coe_plot <- plotCounts(dds, gene = which.min(control_vs_oe$padj),
                      intgroup = "genotype", returnData = TRUE)
koe_plot <- plotCounts(dds, gene = which.min(ko_vs_oe$padj),
                      intgroup = "genotype", returnData = TRUE)


ggplot(cko_plot, aes(x = genotype, y = count, color = genotype)) +
  scale_colour_manual(values = colours)+
  geom_point() +
  theme_minimal() +
  geom_point(shape = 19, size = 2) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))  +
  ggtitle("Control vs KO") +
  scale_y_log10()

ggplot(coe_plot, aes(x = genotype, y = count, color = genotype)) +
  scale_colour_manual(values = colours) +
  geom_point() +
  theme_minimal() +
  geom_point(shape = 19, size = 2) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))  +
  ggtitle("Control vs OE") +
  scale_y_log10()

ggplot(koe_plot, aes(x = genotype, y = count, color = genotype)) +
  scale_colour_manual(values = colours) +
  geom_point() +
  theme_minimal() +
  geom_point(shape = 19, size = 2) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("KO vs OE") +
  scale_y_log10()