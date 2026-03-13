############load libraries
library(dplyr)
library(readr)
library(tidyverse)
##############Editing rate functions
EditingRate_function <-function(sample_list) {
  sample_list$Editing.rate <- NA  
  for (n in 1:nrow(sample_list)) {
    if (sample_list$REF[n] == "A") {
      sample_list$Editing.rate[n] <- sample_list$G[n] / sample_list$Total[n]
    } else {
      sample_list$Editing.rate[n] <- sample_list$C[n] / sample_list$Total[n]
    }
  }
  
  return(sample_list)
}
###filter editing rates
filter_editingrate <- function(df) {
  df %>%
    filter(Editing.rate <= 0.99 & Editing.rate >= 0.01) %>%
    filter(Editing.rate <= 0.49 | Editing.rate >= 0.51)
}
