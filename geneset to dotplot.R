#Geneset to dotplot
#Author: Dr Chris J Smith

#load packages
library(tidyverse)

#load data
control_vs_ko <- read.csv('skin_control_vs_ko_deseq2.csv')
control_vs_oe <- read.csv('skin_control_vs_oe_deseq2.csv')
ko_vs_oe <- read.csv('skin_ko_vs_oe_deseq2.csv')
geneset <- read.csv('geneset.csv')

#filter for geneset and combine comparisons
geneset_control_vs_ko <- control_vs_ko %>%
  inner_join(geneset, by = "Gene") %>%
  mutate(fc =log2FoldChange * -1) %>%
  select(Gene, fc, padj) %>%
  add_column(comparison = "Control vs KO") 

geneset_control_vs_oe <- control_vs_oe %>%
  inner_join(geneset, by = "Gene") %>%
  mutate(fc =log2FoldChange * -1) %>%
  select(Gene, fc, padj) %>%
  add_column(comparison = "Control vs OE")

geneset_ko_vs_oe <- ko_vs_oe %>%
  inner_join(geneset, by = "Gene") %>%
  mutate(fc =log2FoldChange * 1) %>%
  select(Gene, fc, padj) %>%
  add_column(comparison = "KO vs OE")

combined_genelist <- geneset_control_vs_ko %>%
  rbind(geneset_control_vs_oe) %>%
  rbind(geneset_ko_vs_oe) %>%
  arrange(Gene)%>%
  na.omit() 
 
#plot dotplot
ggplot(combined_genelist, aes(comparison, Gene)) +
         geom_point(aes(size = padj,  colour = as.numeric(fc), shape = padj)) + 
  scale_colour_gradient2("log2Fold Change", 
                        breaks = c(-3,-2, -1, 0, 1, 2, 3), 
                        labels = format(c("-3", "-2", "-1","0","1","2","3")),
                        low = "blue", mid = "white", high = "red") +
  scale_shape_binned(name = "Significance Threshold",
                     breaks = c(0.05, 1)) +
  scale_size("padj", trans = 'reverse',
             breaks = c(0.05, 0.1, 0.4), 
             labels = format(c("0.05", "0.10", "0.40")),
             range = c(1, 6)) +
   scale_y_discrete(limits = rev)+
   theme_minimal() +
   theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Geneset")



                        

