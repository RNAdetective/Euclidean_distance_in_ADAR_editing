############load libraries
library(dplyr)
library(readr)
library(stringr)
library(tidyverse)
library(data.table)
############################Format data to calculate pairwise distance 

HC_matrix_data <- data.frame(Coordinate=HC_BL$Coordinate, Editing.rate=HC_BL$Editing.rate,
                             Sample=HC_BL$Run, Gene= HC_BL$snpEff_geneName)

PRO_matrix_data <- data.frame(Coordinate=PRO_BL$Coordinate, Editing.rate=PRO_BL$Editing.rate,
                              Sample=PRO_BL$Run, Gene= PRO_BL$snpEff_geneName)



############################
HC_matrix_data$G.Location <- paste0(HC_matrix_data$Coordinate, ":",HC_matrix_data$Gene)
PRO_matrix_data$G.Location <- paste0(PRO_matrix_data$Coordinate, ":",PRO_matrix_data$Gene)

Count_matrix <- rbind(HC_matrix_data, PRO_matrix_data)
write.csv(Count_matrix, "Editing_Rates_Count_matrix")

Count_matrix$Coordinate <- NULL
Count_matrix$Gene <- NULL

Count_matrix<-   Count_matrix %>%
  pivot_wider(names_from = Sample, values_from = Editing.rate)

Count_matrix<- Count_matrix%>% replace (is.na(.),0)
cm<- as.data.table(Count_matrix)
cm<-as.matrix(cm, rownames = "G.Location")
dim(cm)

############calculate distance
d <- apply(cm, 1, function(x) {
  HC_vals <- x[1:9]
  PRO_vals  <- x[10:27]
  sqrt(sum((HC_vals - PRO_vals)^2))
})
distance_data <- as.data.frame(d)
write.csv(distance_data, "Pairwise_distance.csv")
################################
dGreaterThan1 <- filter(Pairwise_distance, d >= 1,)
PD_risk_GeneList <- left_join(dGreaterThan1, DisGeNET_PDriskGeneList, by="Gene")
PD_risk_GeneList<- PD_risk_GeneList[!is.na(PD_risk_GeneList$Disease_Association_score), ]
write.csv(PD_risk_GeneList, "PD_risk_GeneList.csv")
