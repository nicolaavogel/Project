library(here)
library(vegan)
library(readr)
library(dplyr)
library(caret)
library(stringr)

## Data 

# Checking for overrepresentation of one of the classes. 
is.Forest <- nrow(Data[Data$Forest == 1, ])
is.Forest 
is.not.Forest <- nrow(Data[Data$Forest == 0, ])
is.not.Forest

# 52 to 78

# Model 0: Random Forest on NMDS ordination axes of Bray-Curtis dissimilarities of  transformed "raw" read counts.
# Model 0: NMDS with 6 axis -> used for model building 3 axis (equivalent to the publication.)

fun_ord0 <- metaMDS(t(fun_otu_tab[,1:130]), k=6, try=200, trymax = 1000)
col_dat0 <- data.frame(fun_ord0$points[,1:3])

fungi_0 <- cbind(col_dat0, Data$Forest)
names(fungi_0)[4] <- "Forest"
fungi_0$Forest <- as.factor(fungi_0$Forest)

form_F <- as.formula(paste("Forest ~ .")) 

# Model recreation of the publication -- predictors taken from "DataAnalysis for Publication.R"
# allocated with the stepAIC method and safed in the ResultsBi table 
ctrl <- trainControl(method = "LOOCV", number = 10)
model_pub <- train(form_F, data = fungi_0, method = "qda", trControl = ctrl, metric = "Accuracy")
model_pub
cm <-confusionMatrix(model_pub$pred$pred, model_pub$pred$obs)

all.models <- function(data, form){ 
  
  # All resampling with LOOCV 
  ctrl <- trainControl(method = "LOOCV", number = 10)
  
  # QDA model
  model0 <- train(form, data = data, method = "qda", trControl = ctrl, metric = "Kappa")
  con0 <- confusionMatrix(model0$pred$pred, model0$pred$obs)
  
  # Random Forest model
  model0_RF <- train(form, data = data, method = "rf", trControl = ctrl, metric = "Kappa")
  
  # Only getting the prediction from the best value for Kappa 
  best.mtry <- model0_RF$bestTune[[1]]
  best.model <- model0_RF$pred[model0_RF$pred$mtry == best.mtry, ]

  con0_RF <- confusionMatrix(best.model$pred, best.model$obs)
  
  # K-nearest neighbor model 
  model0_knn <- train(form , data = data, method = "knn", trControl = ctrl, metric = "Kappa")
  
  # Only getting the prediction from the best value for Kappa 
  best.k <- model0_knn$bestTune[[1]]
  best.knn <- model0_knn$pred[model0_knn$pred$k == best.k, ]
  con0_knn <- confusionMatrix(best.knn$pred, best.knn$obs)
  
  # Adding the models and their confusion matrices to a list 
  qda <- list(model0, con0)
  rf <- list(model0_RF, con0_RF)
  knn <- list(model0_knn, con0_knn)
  models.list <- list(qda, rf, knn)

  return(models.list)

}



models_0 <- all.models(fungi_0, form_F)
models_0


# Model 1: Random Forest on NMDS ordination axes of Bray-Curtis dissimilarities of transformed "raw" read counts.
# Model 1: Higher number of predictors (4 instead of 3) taken for model prediction 

fun_ord1 <- metaMDS(t(fun_otu_tab[,1:130]), k=6, try=200, trymax = 2000)
col_dat1 <- data.frame(fun_ord1$points[,1:4])

fungi_1 <- cbind(col_dat1, Data$Forest)
names(fungi_1)[5] <- "Forest"
fungi_1$Forest <- as.factor(fungi_1$Forest)

models_1 <- all.models(fungi_1, form_F)
models_1

mds.data <- data.frame(Sample = rownames(fungi_1), 
                       color = fungi_1$Forest,
                       X = fungi_1[,1],
                       Y = fungi_1[,2])
MDS_Deciduous <- ggplot()+
  geom_text(data = mds.data, aes(x = X, y = Y, label = Sample, color = color)) + 
  theme_bw() 

### Normalisation with sum of row devision and taking the square-root of 4. 
# Model 2 (No scalingof the data; Random Forest on the whole data set)

fun_t <- t(fun_otu_tab[,1:130])
rel_abund_fun_t <- sweep(fun_t, 1, rowSums(fun_t), "/")
rel_abund_fun_t <- rel_abund_fun_t^0.25
fungi_2 <- cbind(rel_abund_fun_t, Data$Forest)
fungi_2 <- as.data.frame(fungi_2)
names(fungi_2)[10491] <- "Forest"
fungi_2$Forest <- as.factor(fungi_2$Forest)


# Adjusting the function for the whole data set: QDA is not suitable for the size of the dataset
two.models <- function(data, form){ 
 
 ctrl <- trainControl(method = "LOOCV", number = 10)
 
 model0_RF <- train(form , data = data, method = "rf", trControl = ctrl, metric = "Kappa")
 
 best.mtry <- model0_RF$bestTune[[1]]
 best.model <- model0_RF$pred[model0_RF$pred$mtry == best.mtry, ]
 
 con0_RF <- confusionMatrix(best.model$pred, best.model$obs)
 
 model0_knn <- train(form , data = data, method = "knn", trControl = ctrl, metric = "Kappa")
 best.k <- model0_knn$bestTune[[1]]
 best.knn <- model0_knn$pred[model0_knn$pred$k == best.k, ]
 con0_knn <- confusionMatrix(best.knn$pred, best.knn$obs)
 

 rf <- list(model0_RF, con0_RF)
 knn <- list(model0_knn, con0_knn)
 models.list <- list(rf, knn)
 
 return(models.list)
 
}


models_2 <- two.models(fungi_2, form_F)
models_2


#### Collapse on genus level // prediction of Forest data 

fungi <- cbind(fun_otu_tab[,1:130], fun_otu_tab$genus)
names(fungi)[131] <- "genus"


fungi_colla <- aggregate( . ~ genus, fungi, sum) 
fungi_colla <- fungi_colla[-990,]    # Removing row 990 - row includes unidentified species and could highly bias the data.


##### Model on collapsed data genus level 
### Normalization method as is model 2.  

colla_test <- t(fungi_colla[,2:131])
fun_t <- colla_test
rel_abund_fun_t <- sweep(fun_t, 1, rowSums(fun_t), "/")
rel_abund_fun_t <- rel_abund_fun_t^0.25
fungi_3 <- cbind(rel_abund_fun_t, Data$Forest)
fungi_3 <- as.data.frame(fungi_3)
names(fungi_3)[1038] <- "Forest"
fungi_3$Forest <- as.factor(fungi_3$Forest)

models_3 <- two.models(fungi_3, form_F)
models_3


# Using the collapsed data with the same scaling normalization as in model 1:

col_fun_sca <- metaMDS(t(fungi_colla[,2:131]), k = 4, try = 200, trymax = 1000)
col_dat_col <- data.frame(col_fun_sca$points[,1:4])
fungi_4 <- cbind(col_dat_col, Data$Forest)
names(fungi_4)[5] <- "Forest"
fungi_4$Forest <- as.factor(fungi_4$Forest)

models_4 <- all.models(fungi_4, form_F)
models_4

# Adding the plant dataset as predictor variables with its' own NMDS ordination of Bray-Curtis dissimilarity. Fungi is used as in model 1. 

pla_otu_tab <- read.table(here::here("data","plant_table_final.txt"), sep="\t", header = T)
pla_ord <- metaMDS(t(euk_otu_tab[,1:130]), k=4, try= 200, trymax = 1000)
col_dat_pla <- data.frame(pla_ord$points[,1:4])


fungi_5 <- cbind(col_dat1, col_dat_pla, Data$Forest)
names(fungi_5) <- c("fung1","fung2","fung3","fung4", "plant1", "plant2", "plant3", "plant4")
names(fungi_5)[9] <- "Forest"
fungi_5$Forest <- as.factor(fungi_5$Forest)

models_5 <- all.models(fungi_5,form_F)
models_5



#############################################################
# Predicting different habitats: Atlantic/Oak/Willow 
library(DMwR)

#### ATLANTIC 
# Checking for overrepresentation of one of the classes. 
is.Atlantic <- nrow(Data[Data$Atlantic == 1, ])
is.Atlantic 
is.not.Atlantic <- nrow(Data[Data$Atlantic == 0, ])
is.not.Atlantic

# 35 to 95 

fungi_5_A <- cbind(col_dat1, col_dat_pla, Data$Atlantic)
names(fungi_5_A) <- c("fung1","fung2","fung3","fung4", "plant1", "plant2", "plant3", "plant4")
names(fungi_5_A)[9] <- "Atlantic"
fungi_5_A$Atlantic <- as.factor(fungi_5_A$Atlantic)

form_A <- as.formula(paste("Atlantic ~ .")) 


models_5_A <- all.models(fungi_5_A, form_A)
models_5_A


# Data adjustment with SMOTE method 
newData <- SMOTE(Atlantic ~. , data = fungi_5_A, perc.over = 300, perc.under = 150)
table(newData$Atlantic)


newData$Atlantic <- as.factor(newData$Atlantic)

models_1_A_new <- all.models(newData, form_A)
models_1_A_new

#### OAK
# Checking for overrepresentation of one of the classes. 
is.Oak <- nrow(Data[Data$Oak == 1, ])
is.Oak
is.not.Oak <- nrow(Data[Data$Oak == 0, ])
is.not.Oak

# 11 to 119 

fungi_5_O <- cbind(col_dat1, col_dat_pla, Data$Oak)
names(fungi_5_O) <- c("fung1","fung2","fung3","fung4", "plant1", "plant2", "plant3", "plant4")
names(fungi_5_O)[9] <- "Oak"
fungi_5_O$Oak <- as.factor(fungi_5_O$Oak)

form_O <- as.formula(paste("Oak ~ .")) 


models_5_O <- all.models(fungi_5_O, form_O)
models_5_O


#######################

ggplot(data = fungi_5_O, aes(fung2, plant4, color = Oak)) + 
 geom_point() + 
 theme_bw()


ggplot(data = newData_O, aes(fung2, plant4, color = Oak)) + 
 geom_point() +
 theme_bw() 


# Data adjustment with SMOTE method 
newData_O <- SMOTE(Oak ~. , data = fungi_5_O, perc.over = 500, perc.under = 150)
table(newData_O$Oak)


newData_O$Oak <- as.factor(newData_O$Oak)

models_5_O_new <- all.models(newData_O, form_O)
models_5_O_new

varImp(models_5[[3]][[1]])


#### WILLOW 
# Checking for overrepresentation of one of the classes. 
is.Willow <- nrow(Data[Data$Willow == 1, ])
is.Willow
is.not.Willow <- nrow(Data[Data$Willow == 0, ])
is.not.Willow

# 10 to 120

fungi_5_W <- cbind(col_dat1,col_dat_pla, Data$Willow)
names(fungi_5_W) <- c("fung1","fung2","fung3","fung4", "plant1", "plant2", "plant3", "plant4")
names(fungi_5_W)[9] <- "Willow"
fungi_5_W$Willow <- as.factor(fungi_5_W$Willow)

form_W <- as.formula(paste("Willow ~ .")) 

models_5_W <- all.models(fungi_5_W, form_W)
models_5_W

# Data adjustment with SMOTE method 
newData_W <- SMOTE(Willow ~. , data = fungi_5_W, perc.over = 500, perc.under = 150)
table(newData_W$Willow)


newData_W$Willow <- as.factor(newData_W$Willow)

models_5_W_new <- all.models(newData_W, form_W)
models_5_W_new
