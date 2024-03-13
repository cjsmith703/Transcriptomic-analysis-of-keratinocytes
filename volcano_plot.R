#Volcano plot from RNAseq
#Author: Dr Chris J Smith

#load packages
library(tidyverse)

#read deseq2 file
deseq2 <- read.csv('ko_vs_oe_deseq2.csv')

#create volcano dataframe
volcano <- deseq2 %>%
  arrange(padj) %>%
  drop_na() %>%
  select(Gene, log2FoldChange, pvalue, padj) %>%
  mutate(logp = -log10(padj))

#set up/down regulated
volcano$diffexpressed <- "Not Sig"
volcano$diffexpressed[volcano$log2FoldChange > 1 & volcano$padj < 0.05] <- "Up Regulated"
volcano$diffexpressed[volcano$log2FoldChange < -1 & volcano$padj < 0.05] <- "Down Regulated"

#set colours
mycolours <- c("#FF0000", "#0000FF", "grey")

names(mycolours) <- c("Up Regulated", "Down Regulated", "Not Sig")

#plot volcano plot with top 5 genes labelled
ggplot(volcano, aes(log2FoldChange, logp, col = diffexpressed)) +
   theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "right") +
  scale_colour_manual(values = mycolours) +
  xlim(-6, 6) +
  xlab("Log2 Fold Change") +
  ylab("-Log10(P adj)") +
  geom_text(data = head(volcano, 5), aes(label = Gene, colour = "black"), nudge_x = 1, show.legend = FALSE ) +
  geom_jitter(size=1)

