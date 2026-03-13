####################Random Forest##############
library(dplyr)     
library(ggplot2)   
library(rsample)  
library(broom)
library(randomForest)
library(tidyverse)
library(caret)
########################prepare test and train data

data_test <- HC_prodromalPD_test_data
data_train <- HC_prodromalPD_training_data

data_train<-   data_train %>%
  pivot_wider(names_from = G.Location, values_from = Editing.rate)

data_train <- data_train%>% replace (is.na(.),0)

data_train$Condition <- as.factor(data_train$Condition)


data_train$Sample <- NULL

colnames(data_train)

colnames(data_train) <- gsub("-", "", colnames(data_train))
#########data test

data_test<-   data_test %>%
  pivot_wider(names_from = G.Location, values_from = Editing.rate)




data_test <- data_test%>% replace (is.na(.),0)
data_test$Condition <- as.factor(data_test$Condition)

colnames(data_test) <- gsub("-", "", colnames(data_test))



##############################

common_colnames <- intersect(colnames(data_train), colnames(data_test))


# Filter both data frames to keep only common columns
data_train <- data_train[, common_colnames]
data_test <- data_test[, common_colnames]


#######################RF model
set.seed(42)

model <- randomForest(as.factor(Condition) ~., data = data_train, ntree=500)
model

Actuall_levels <- data_test


data_test$Condition <- NULL


###if no TYpe is selected default is "response"
prediction <- predict(model, newdata = data_test, type="response")
prediction

conf_matrix <- confusionMatrix(prediction, Actuall_levels$Condition)
conf_matrix$table

######################################Prodromal-PD and PD 
data_train2 <- PD_prodromal_PD_training_data

data_test2 <- PD_prodromal_PD_test_data

data_train2<-   data_train2 %>%
  pivot_wider(names_from = G.Location, values_from = Editing.rate)


data_train2 <- data_train2%>% replace (is.na(.),0)

data_train2$Condition <- as.factor(data_train2$Condition)
colnames(data_train2)<- gsub("-","", colnames(data_train2) )

#########data test

data_test2<-   data_test2 %>%
  pivot_wider(names_from = G.Location, values_from = Editing.rate)


data_test2 <- data_test2%>% replace (is.na(.),0)

data_test2$Condition <- as.factor(data_test2$Condition)

common_colnames2 <- intersect(colnames(data_train2), colnames(data_test2))


# Filter both data frames to keep only common columns
data_train2 <- data_train2[, common_colnames2]
data_test2 <- data_test2[, common_colnames2]
data_train2$Sample <- NULL

###############RF model
set.seed(42)
model2 <- randomForest(as.factor(Condition) ~., data = data_train2, ntree=500)
model2
data_test2$Sample <- NULL


Actuall_levels2 <- data_test2


data_test2$Condition <- NULL



prediction2 <- predict(model2, newdata = data_test2, type="response")
prediction2

conf_matrix2 <- confusionMatrix(prediction2, Actuall_levels2$Condition)

########################Notes##############
###Model's performance may vary with iterations as implied by large CI prediction intervals 

