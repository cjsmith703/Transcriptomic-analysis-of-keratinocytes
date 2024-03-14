#GSEA to Dotplot 
#Author: Dr Chris J Smith

#load packages
library(tidyverse)

#load data
control_vs_ko <- read.csv("data/GSEA_control_vs_ko.csv")
control_vs_oe <- read.csv("data/GSEA_control_vs_oe.csv")
ko_vs_oe <- read.csv("data/GSEA_ko_vs_oe.csv")

#filter, slice top 10 genesets and combine data
top_control_vs_ko <- control_vs_ko %>%
  filter(FDR < 0.05) %>%
  arrange(desc(ABS_NES)) %>%
  mutate(new_NES = NES *-1) %>%
  slice(1:10) %>%
  add_column(comparison = "Control vs KO")

top_control_vs_oe <- control_vs_oe %>%
  filter(FDR < 0.05) %>%
  arrange(desc(ABS_NES)) %>%
  mutate(new_NES = NES *-1) %>%
  slice(1:10) %>%
  add_column(comparison = "Control vs OE")

top_ko_vs_oe <- ko_vs_oe%>%
  filter(FDR < 0.05) %>%
  arrange(desc(ABS_NES)) %>%
  mutate(new_NES = NES *-1) %>%
  slice(1:10) %>%
  add_column(comparison = "KO vs OE")

gsea <- top_control_vs_ko %>%
  rbind(top_control_vs_oe) %>%
  rbind(top_ko_vs_oe) %>%
  mutate(GS = gsub("KEGG_", "", NAME))%>%
  mutate(GS = gsub("_", " ", GS))

#plot dotplot
ggplot(gsea, aes(comparison, GS))+
  geom_point(aes(colour = as.numeric(new_NES), size = FDR))+
  theme_minimal()+
  scale_size("FDR", trans = 'reverse')+
  scale_colour_gradient2("NES", 
                         breaks = c(-3,-2, -1, 0, 1, 2, 3), 
                         labels = format(c("-3", "-2", "-1","0","1","2","3")),
                         low = "blue", mid = "white", high = "red") +
  scale_y_discrete(limits = rev)+
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("GSEA")


