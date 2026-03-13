library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)
library(ggsci)
library(ggstatsplot)
library(ggsignif)
library(tidyr)
###########Density plots
HC_BL <- HC_BL[grepl("protein_coding", HC_BL$snpEff_Transcript_BioType, fixed = TRUE), ]
HC_density_plot <- ggplot(HC_BL, aes(x=Editing.rate))+
  geom_density(color="black", fill="lightblue") + theme_classic()


HC_density_plot <- HC_density_plot + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.5)))
HC_density_plot <- HC_density_plot+theme(axis.text.y = element_text(face="bold"), axis.title.y = element_text(face="bold"),axis.title.x = element_text(face="bold"),
                                         axis.text.x = element_text(face="bold"))


HC_density_plot




PRO_BL <- PRO_BL[grepl("protein_coding", PRO_BL$snpEff_Transcript_BioType, fixed = TRUE), ]
PRO_density_plot <- ggplot(PRO_BL, aes(x=Editing.rate))+
  geom_density(color="black", fill="purple1") + theme_classic()

PRO_density_plot <- PRO_density_plot + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.5)))
PRO_density_plot <- PRO_density_plot+theme(axis.text.y = element_text(face="bold"), axis.title.y = element_text(face="bold"),axis.title.x = element_text(face="bold"),
                                           axis.text.x = element_text(face="bold"))


PRO_density_plot

Density_plot <- plot_grid(
  HC_density_plot, PRO_density_plot,nrow = 1, labels = c("A", "B") )

Density_plot  


ggsave(
  filename = "DENSITY_PLOT.png",
  plot = Density_plot,
  width = 10,
  height = 4,
  units = "in",
  dpi = 600
)
####################ADAR_gene_plot
names(normalised_counts_gene)[1] <- "Gene"
normalised_counts_gene <- normalised_counts_gene %>% separate(col = Gene,into = c("Gene_ID","Gene") , 
                                                              sep="\\|")

ADAR_gene_Counts <-  normalised_counts_gene[grepl("ADAR", normalised_counts_gene$Gene, fixed = TRUE), ]
ADAR1_counts <- ADAR_gene_Counts[grepl("ENSG00000160710", ADAR_gene_Counts$Gene_ID, fixed = TRUE),]
ADAR2_counts <- ADAR_gene_Counts[grepl("ENSG00000185736", ADAR_gene_Counts$Gene_ID, fixed = TRUE),]
ADAR3_counts <- ADAR_gene_Counts[grepl("ENSG00000185736", ADAR_gene_Counts$Gene_ID, fixed = TRUE),]
names(ADAR1_counts)[3:11] <- "Healthy_Controls"
names(ADAR1_counts)[12:29] <- "Prodromal_PD"
ADAR1_counts$Gene_ID <- NULL
ADAR1_counts$Gene <- NULL
ADAR1_counts <- ADAR1_counts %>%
  pivot_longer(
    cols = everything(),        
    names_to = "Condition",
    values_to = "Gene Counts"
  )

p1 <- ggbetweenstats(
  data = ADAR1_counts,
  x = Condition,
  y = `Gene Counts`,
  type = "nonparametric",
  results.subtitle = TRUE,
  pairwise.display = "all",
  pairwise.comparison= TRUE,
  p.adjust.method = "BH",
  , pairwise.annotation = "p.value", 
  conf.level = 0.95,
  title = "ADAR1",
  package = "ggsci",
  palette = "nrc_npg") + theme_classic()

p1


names(ADAR2_counts)[3:11] <- "Healthy_Controls"
names(ADAR2_counts)[12:29] <- "Prodromal_PD"
ADAR2_counts$Gene_ID <- NULL
ADAR2_counts$Gene <- NULL
ADAR2_counts <- ADAR2_counts %>%
  pivot_longer(
    cols = everything(),        
    names_to = "Condition",
    values_to = "Gene Counts"
  )

p2 <- ggbetweenstats(
  data = ADAR2_counts,
  x = Condition,
  y = `Gene Counts`,
  type = "nonparametric",
  results.subtitle = TRUE,
  pairwise.display = "all",
  pairwise.comparison= TRUE,
  p.adjust.method = "BH",
  , pairwise.annotation = "p.value", 
  conf.level = 0.95,
  title = "ADAR2",
  package = "ggsci",
  palette = "nrc_npg") + theme_classic()

p2


names(ADAR3_counts)[3:11] <- "Healthy_Controls"
names(ADAR3_counts)[12:29] <- "Prodromal_PD"
ADAR3_counts$Gene_ID <- NULL
ADAR3_counts$Gene <- NULL
ADAR3_counts <- ADAR3_counts %>%
  pivot_longer(
    cols = everything(),        
    names_to = "Condition",
    values_to = "Gene Counts"
  )

p3 <- ggbetweenstats(
  data = ADAR3_counts,
  x = Condition,
  y = `Gene Counts`,
  type = "nonparametric",
  results.subtitle = TRUE,
  pairwise.display = "all",
  pairwise.comparison= TRUE,
  p.adjust.method = "BH",
  , pairwise.annotation = "p.value", 
  conf.level = 0.95,
  title = "ADAR3",
  package = "ggsci",
  palette = "nrc_npg") + theme_classic()

p3

ADAR_gene_plot <- combine_plots(
  list(p1, p2,p3),
  plotgrid.args = list(nrow = 1)
  
  
)

ADAR_gene_plot

ggsave("ADAR_gene_plot.svg",plot=ADAR_gene_plot, device=svglite::svglite, width=8000, height=5000, scale=0.8, units="px", bg="white")
###################################
######################Confusion matrix plots
#####################################
CM1 <- as.data.frame(conf_matrix$table)
metrics_text <- paste0(
  "Classification Report\n",
  "-------------------------\n",
  "Accuracy: 0.7778\n",
  "95% CI: (0.3999, 0.9719)\n",
  "Sensitivity: 0.3333\n",
  "Specificity:1.0000\n",
  "Balanced Accuracy: 0.6667\n\n",
  "Positive Class: Healthy Control"
)

plt1 <- ggplot(CM1, aes(Prediction, Reference, fill = Freq)) +
  
  geom_tile(color = "white", linewidth = 1.2) +
  
  geom_text(aes(label = Freq),
            size = 7,
            fontface = "bold") +
  
  scale_fill_gradient(low = "#c6dbef", high = "#08306b",
                      name = "Frequency") +
  
  labs(
    x = "True Class",
    y = "Predicted Class",
    title = "Confusion Matrix"
  ) +
  
  scale_x_discrete(labels = c( "Prodromal PD","Healthy Control")) +
  scale_y_discrete(labels = c("Healthy Control","Prodromal PD" )) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.line = element_blank()
  )

plt1

classification_report <- ggdraw() +
  draw_label(metrics_text,
             x = 0,
             hjust = 0,
             fontface = "bold",
             size = 13) +
  theme(
    plot.background = element_rect(fill = "white", color = NA
    )
  )

CM_HC_PRO_PD <- plot_grid(
  plt1,
  classification_report,
  ncol = 2,
  rel_widths = c(2, 1)
)


###################################
CM2 <- as.data.frame(conf_matrix2$table)
metrics_text2 <- paste0(
  "Classification Report\n",
  "-------------------------\n",
  "Accuracy : 0.8333\n",
  "95% CI: (0.5159, 0.9791)\n",
  "Sensitivity:0.8333\n",
  "Specificity:  0.8333\n",
  "Balanced Accuracy: 0.8333\n\n",
  "Positive Class: PD"
)

PLOT2 <- ggplot(CM2, aes(Prediction, Reference, fill = Freq)) +
  
  geom_tile(color = "white", linewidth = 1.2) +
  
  geom_text(aes(label = Freq),
            size = 7,
            fontface = "bold") +
  
  scale_fill_gradient(low = "lightpink", high = "deeppink3",
                      name = "Frequency") +
  
  labs(
    x = "True Class",
    y = "Predicted Class",
    title = "Confusion Matrix"
  ) +
  
  scale_x_discrete(labels = c( "Prodromal PD","PD")) +
  scale_y_discrete(labels = c("PD","Prodromal PD" )) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.line = element_blank()
  )

PLOT2

classification_report2 <- ggdraw() +
  draw_label(metrics_text2,
             x = 0,
             hjust = 0,
             fontface = "bold",
             size = 13) +
  theme(
    plot.background = element_rect(fill = "white", color = NA
    )
  )

CM_PD_PRO_PD <- plot_grid(
  PLOT2,
  classification_report2,
  ncol = 2,
  rel_widths = c(2, 1)
)

CM_plots <- plot_grid(CM_HC_PRO_PD,CM_PD_PRO_PD ,nrow = 1, labels = c("A", "B") )



ggsave(
  filename = "CM_plots.png",
  plot = CM_plots,
  width = 17,
  height = 6,
  units = "in",
  dpi = 600
)

##################GOanalysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(GOSemSim)
library(readr)
library(enrichplot)
library(DOSE)
library(ggrepel)
library(gtable)


ego_all <- enrichGO(gene =shared_and_uniqueBiomarker_gene_list$All_genes,
                    OrgDb = org.Hs.eg.db,
                    keyType = 'SYMBOL',
                    ont = "ALL",
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05, 
                    maxGSSize   = 50,
                    readable = TRUE) 


ego_all_genes <- ego_all@result
write.csv(ego_all_genes, "ego_all_genes.csv")

p_facet <- dotplot(
  ego_all,
  x = "Count",
  showCategory = 15,
  split = "ONTOLOGY"
) +
   enrichplot::autofacet(by = "row", scales = "free")

all_gene_plot <- p_facet + scale_y_discrete()+theme(axis.text.y = element_text(face="bold"), axis.text.x = element_text(face="bold"),legend.text = element_text(face = "bold"),
                                                    legend.title = element_text(face = "bold")) + labs(x="Gene count")


ego_PRO <- enrichGO(gene =shared_and_uniqueBiomarker_gene_list$Unique_to_PRO,
                    OrgDb = org.Hs.eg.db,
                    keyType = 'SYMBOL',
                    ont = "ALL",
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05, 
                    maxGSSize   = 50,
                    readable = TRUE) 


ego_Pro_unique_results <- ego_PRO@result
write.csv(ego_Pro_unique_results, "ego__unique_pro_results.csv")

p_facet2 <- dotplot(
  ego_PRO,
  x = "Count",
  showCategory = 15,
  split = "ONTOLOGY"
) +
  enrichplot::autofacet(by = "row", scales = "free")

PRO_gene_plot <- p_facet2 + scale_y_discrete()+theme(axis.text.y = element_text(face="bold"), axis.text.x = element_text(face="bold"),legend.text = element_text(face = "bold"),
                                                    legend.title = element_text(face = "bold")) + labs(x="Gene count")






GO_plot <- plot_grid( all_gene_plot,PRO_gene_plot,nrow = 1, labels = c("A", "B") )

GO_plot

ggsave(
  filename = "GO_plot.png",
  plot =GO_plot,
  width = 22,
  height = 10,
  units = "in",
  dpi = 800
)











