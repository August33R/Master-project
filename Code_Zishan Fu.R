rm(list = ls())
setwd("E:/you know files/IC/PD project/Data2")

library(devtools)
library(dplyr)
library(ape)
library(phytools)
library(mallorn)  ## PD calculation based on numeric prob of extinction
library(pastecs)  ## Package for Analysis of Space-Time Ecological Series (stat.desc Descriptive statistics on a data frame or time series)
library(stats)
install.packages('psych')
library(psych)
library(ggplot2)

## README - the codes are divided into sections. Please use the small triangles to pull in and out.

## preperation---------------------------------------------------------
## Probability transformation----------------------

## Reading and editing trees
Coral_tree <- read.nexus(file = "Reef_randChrono.tre")

tips_of_syn <- c("ACR_Acropora_muricata", "DEN_Turbinaria_crater", "FAV_Colpophyllia_amaranthus",
                 "FAV_Colpophyllia_breviserialis", "MEA_Meandrina_jacksoni")
Coral_tree_nosyn <- lapply(Coral_tree,drop.tip,tip=tips_of_syn) ##Coral tree without synonyms
class(Coral_tree_nosyn)<-"multiPhylo"
writeNexus(Coral_tree_nosyn, file="E:/you know files/IC/PD project/Data2/Coral_tree_nosyn2.tre")

tips_of_noRL <- c('ACR_Acropora_intermedia','ACR_Astreopora_acroporina','ACR_Astreopora_cenderawasih',
                  'ACR_Astreopora_montiporina','ACR_Isopora_meridiana','FAV_Favia_gravida',
                  'FAV_Goniastrea_thecata','FUN_Lithophyllon_puishani','FUN_Pleuractis_gravis',
                  'FUN_Podabacia_kunzmanni','POR_Porites_randalli','SID_Craterastrea_levis')
Coral_tree_mRL <- lapply(Coral_tree_nosyn, drop.tip,tip=tips_of_noRL) ##Coral tree without the species not in the RL
class(Coral_tree_mRL)<-"multiPhylo"
writeNexus(Coral_tree_mRL, file="E:/you know files/IC/PD project/Data2/Coral_tree_mRL2.tre")

## Converting Endangered Species Categories to Probabilities of Extinction----------------
convertProb <- function(Coral_RL){
  Coral_RL %>%
    mutate(Current_Issac = recode(Coral_RL$Current_Cat, 'LC' = 0.025, 'NT' = 0.05, 'VU' = 0.1, 'EN' = 0.2, 'CR' = 0.4,'EX' = 1, 'EW' = 1))%>%
    mutate(Current_IUCN40 = recode(Coral_RL$Current_Cat, 'LC' = 0.00004, 'NT' = 0.00401, 'VU' = 0.04127, 'EN' = 0.36000, 'CR' = 0.93750,'EX' = 1, 'EW' = 1)) %>%
    mutate(Current_IUCN50 = recode(Coral_RL$Current_Cat, 'LC' = 0.00005, 'NT' = 0.00501, 'VU' = 0.05132, 'EN' = 0.42757, 'CR' = 0.96875,'EX' = 1, 'EW' = 1)) %>%
    mutate(Current_IUCN90 = recode(Coral_RL$Current_Cat, 'LC' = 0.00009, 'NT' = 0.00900, 'VU' = 0.09047, 'EN' = 0.63364, 'CR' = 0.99805,'EX' = 1, 'EW' = 1))%>%
    mutate(Current_IUCN100 = recode(Coral_RL$Current_Cat, 'LC' = 0.00010, 'NT' = 0.01000, 'VU' = 0.10000, 'EN' = 0.67232, 'CR' = 0.99902,'EX' = 1, 'EW' = 1))%>%
    mutate(Current_IUCN490 = recode(Coral_RL$Current_Cat, 'LC' = 0.00049, 'NT' = 0.04805, 'VU' = 0.40326, 'EN' = 0.99578, 'CR' = 1,'EX' = 1, 'EW' = 1))%>%
    mutate(Current_IUCN500 = recode(Coral_RL$Current_Cat, 'LC' = 0.00050, 'NT' = 0.04901, 'VU' = 0.40951, 'EN' = 0.99622, 'CR' = 1,'EX' = 1, 'EW' = 1)) %>%
    mutate(Current_Pess = recode(Coral_RL$Current_Cat, 'LC' = 0.2, 'NT' = 0.4, 'VU' = 0.8, 'EN' = 0.9, 'CR' = 0.99,'EX' = 1, 'EW' = 1))%>%
    mutate(Current_Binary = recode(Coral_RL$Current_Cat, 'LC' = 0, 'NT' = 0, 'VU' = 1, 'EN' = 1, 'CR' = 1,'EX' = 1, 'EW' = 1))%>%
    
    mutate(Past_Issac = recode(Coral_RL$Past_Cat, 'LC' = 0.025, 'NT' = 0.05, 'VU' = 0.1, 'EN' = 0.2, 'CR' = 0.4,'EX' = 1, 'EW' = 1))%>%
    mutate(Past_IUCN40 = recode(Coral_RL$Past_Cat, 'LC' = 0.00004, 'NT' = 0.00401, 'VU' = 0.04127, 'EN' = 0.36000, 'CR' = 0.93750,'EX' = 1, 'EW' = 1)) %>%
    mutate(Past_IUCN50 = recode(Coral_RL$Past_Cat, 'LC' = 0.00005, 'NT' = 0.00501, 'VU' = 0.05132, 'EN' = 0.42757, 'CR' = 0.96875,'EX' = 1, 'EW' = 1)) %>%
    mutate(Past_IUCN90 = recode(Coral_RL$Past_Cat, 'LC' = 0.00009, 'NT' = 0.00900, 'VU' = 0.09047, 'EN' = 0.63364, 'CR' = 0.99805,'EX' = 1, 'EW' = 1)) %>%
    mutate(Past_IUCN100 = recode(Coral_RL$Past_Cat, 'LC' = 0.00010, 'NT' = 0.01000, 'VU' = 0.10000, 'EN' = 0.67232, 'CR' = 0.99902,'EX' = 1, 'EW' = 1))%>% 
    mutate(Past_IUCN490 = recode(Coral_RL$Past_Cat, 'LC' = 0.00049, 'NT' = 0.04805, 'VU' = 0.40326, 'EN' = 0.99578, 'CR' = 1,'EX' = 1, 'EW' = 1)) %>%
    mutate(Past_IUCN500 = recode(Coral_RL$Past_Cat, 'LC' = 0.00050, 'NT' = 0.04901, 'VU' = 0.40951, 'EN' = 0.99622, 'CR' = 1,'EX' = 1, 'EW' = 1)) %>%
    mutate(Past_Pess = recode(Coral_RL$Past_Cat, 'LC' = 0.2, 'NT' = 0.4, 'VU' = 0.8, 'EN' = 0.9, 'CR' = 0.99,'EX' = 1, 'EW' = 1))%>%
    mutate(Past_Binary = recode(Coral_RL$Past_Cat, 'LC' = 0, 'NT' = 0, 'VU' = 1, 'EN' = 1, 'CR' = 1,'EX' = 1, 'EW' = 1))
}

Coral_Prob_nosyn2 <- convertProb(Coral_RL)

write.csv(Coral_Prob, file = "E:/you know files/IC/PD project/Data2/Coral_Prob_nosyn2.csv")



## data imputation (missing data and DD species)-------------------------

dataimptCat <- function(Coral_RL){
  Coral_RL_DD.index <- which((Coral_RL$Current_Cat=='DD')==T) 
  Coral_RL.nonDD <- Coral_RL[-Coral_RL_DD.index,]  
  Coral_RL_DD <- Coral_RL[Coral_RL_DD.index,]
  Coral_RL_DD$Current_Cat <-sample(Coral_RL.nonDD$Current_Cat,nrow(Coral_RL_DD),replace = T)
  Coral_RL_DD$Past_Cat <-sample(Coral_RL.nonDD$Past_Cat,nrow(Coral_RL_DD),replace = T)
  Coral_RL_edited <- rbind(Coral_RL_DD, Coral_RL.nonDD)
}

dataimpt <- function(Coral_Prob_edited){
  Coral_Prob_na.index <- which(is.na(Coral_Prob_edited$Current_Issac)==T)
  Coral_Prob.clean <- Coral_Prob_edited[-Coral_Prob_na.index,]  
  Coral_Prob.na<-Coral_Prob_edited[Coral_Prob_na.index,]
  
  Coral_Prob.na$Current_Issac <-sample(Coral_Prob.clean$Current_Issac,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Current_IUCN40 <-sample(Coral_Prob.clean$Current_IUCN40,nrow(Coral_Prob.na),replace = T)
  Coral_Prob.na$Current_IUCN50 <-sample(Coral_Prob.clean$Current_IUCN50,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Current_IUCN90 <-sample(Coral_Prob.clean$Current_IUCN90,nrow(Coral_Prob.na),replace = T)
  Coral_Prob.na$Current_IUCN100 <-sample(Coral_Prob.clean$Current_IUCN100,nrow(Coral_Prob.na),replace = T)	
  Coral_Prob.na$Current_IUCN490 <-sample(Coral_Prob.clean$Current_IUCN490,nrow(Coral_Prob.na),replace = T)
  Coral_Prob.na$Current_IUCN500 <-sample(Coral_Prob.clean$Current_IUCN500,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Current_Pess <-sample(Coral_Prob.clean$Current_Pess,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Current_Binary <-sample(Coral_Prob.clean$Current_Binary,nrow(Coral_Prob.na),replace = T)		
  
  Coral_Prob.na$Past_Issac <-sample(Coral_Prob.clean$Past_Issac,nrow(Coral_Prob.na),replace = T)	
  Coral_Prob.na$Past_IUCN40<-sample(Coral_Prob.clean$Past_IUCN40,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_IUCN50<-sample(Coral_Prob.clean$Past_IUCN50,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_IUCN90<-sample(Coral_Prob.clean$Past_IUCN90,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_IUCN100 <-sample(Coral_Prob.clean$Past_IUCN100,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_IUCN490<-sample(Coral_Prob.clean$Past_IUCN490,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_IUCN500 <-sample(Coral_Prob.clean$Past_IUCN500,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_Pess <-sample(Coral_Prob.clean$Past_Pess,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_Binary <-sample(Coral_Prob.clean$Past_Binary,nrow(Coral_Prob.na),replace = T)		
  
  rbind(Coral_Prob.clean, Coral_Prob.na)
}

dataimptBi <- function(Prob_binary_edited){
  Prob_binary_na.index <- which(is.na(Prob_binary_edited$Current_Binary)==T)
  Prob_binary.clean <- Prob_binary_edited[-Prob_binary_na.index,]  
  Prob_binary.na<-Prob_binary_edited[Prob_binary_na.index,]
  
  Prob_binary.na$Current_Binary <-sample(Prob_binary.clean$Current_Binary,nrow(Prob_binary.na),replace = T)		
  Prob_binary.na$Past_Binary <-sample(Prob_binary.clean$Past_Binary,nrow(Prob_binary.na),replace = T)	
  rbind(Prob_binary.clean, Prob_binary.na)
}

## Convert list to numeric matrix----------------------------------
convertFormat <- function(Prob_input){
  Prob_matrix <- select(Prob_input, -"Species",-"Current_Cat",-"Past_Cat",-"Current_Binary",-"Past_Binary")
  Prob_Num <- matrix(data = NA, nrow = dim(Prob_matrix)[1], ncol = dim(Prob_matrix)[2])
  for (i in 1:dim(Prob_matrix)[2]) {Prob_Num[,i] <- c(as.numeric(Prob_matrix[[i]]))}
  colnames(Prob_Num) <- c('Current_Issac','Current_IUCN40','Current_IUCN50','Current_IUCN90','Current_IUCN100','Current_IUCN490','Current_IUCN500','Current_Pess',
                          'Past_Issac','Past_IUCN40','Past_IUCN50','Past_IUCN90','Past_IUCN100','Past_IUCN490','Past_IUCN500','Past_Pess')
  rownames(Prob_Num) <- Prob_input$Species
  return(Prob_Num)
}


## Hedge's g---------------------
Hedgesg <- function(x,y){
  mean1 <- mean(x)
  mean2 <- mean(y)
  sd1 <- sd(x)
  sd2 <- sd(y)
  poolSD <- ((((length(x)-1)*sd1^2)+((length(y)-1)*sd2^2))/(length(x)+length(y)-2))^(1/2)
  g <- (mean1-mean2)/poolSD
}


mean1 <- mean(Coral_Past_Issac)
mean2 <- mean(Random_Past_Issac)
sd1 <- sd(Coral_Past_Issac)
sd2 <- sd(Random_Past_Issac)
poolSD <- ((((length(Coral_Past_Issac)-1)*sd1^2)+((length(Random_Past_Issac)-1)*sd2^2))/(length(Coral_Past_Issac)+length(Random_Past_Issac)-2))^(1/2)
g <- (mean1-mean2)/poolSD

##---------------------------------------------------------------------
## pre-analysis--------------------------
setwd("E:/you know files/IC/PD project/Data2")

describe.by(Coral_Current_Issac)
describe.by(Coral_Current_IUCN500)

Coral_Prob_mRL <- read.csv(file = "E:/you know files/IC/PD project/Data2/Coral_Prob_mRL.csv")

##Coral tree without the species not in the RL
tips_of_noRL <- c('ACR_Acropora_intermedia','ACR_Astreopora_acroporina','ACR_Astreopora_cenderawasih',
                  'ACR_Astreopora_montiporina','ACR_Isopora_meridiana','FAV_Favia_gravida',
                  'FAV_Goniastrea_thecata','FUN_Lithophyllon_puishani','FUN_Pleuractis_gravis',
                  'FUN_Podabacia_kunzmanni','POR_Porites_randalli','SID_Craterastrea_levis')
Coral_tree_nosyn2 <- read.nexus(file = "E:/you know files/IC/PD project/Data2/Coral_tree_nosyn2.tre")
Coral_tree_mRL <- lapply(Coral_tree_nosyn, drop.tip,tip=tips_of_noRL)
class(Coral_tree_mRL)<-"multiPhylo"


## PD calculation for 1000 trees
Coral_tree <- Coral_tree_mRL

Coral_PD.mrl <- numeric(1000)
for (j in 1:1000) {
  Coral_Prob_input <- dataimpt(Coral_Prob_mRL)
  Coral_Prob_Num <- convertFormat(Coral_Prob_input)
  Coral_PD.mrl[j] <- ePD(tree = Coral_tree[[j]], probabilities.tips.present.matrix = Coral_Prob_Num)
}

## Results
## Issac
Coral_Current_Issac.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Current_Issac.mrl[k] <- Coral_PD.mrl[[k]][1,1]}
Coral_Past_Issac.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Past_Issac.mrl[k] <- Coral_PD.mrl[[k]][9,1]}

## IUCN40
Coral_Current_IUCN40.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Current_IUCN40.mrl[k] <- Coral_PD.mrl[[k]][2,1]}
Coral_Past_IUCN40.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Past_IUCN40.mrl[k] <- Coral_PD.mrl[[k]][10,1]}

## IUCN50
Coral_Current_IUCN50.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Current_IUCN50.mrl[k] <- Coral_PD.mrl[[k]][3,1]}
Coral_Past_IUCN50.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Past_IUCN50.mrl[k] <- Coral_PD.mrl[[k]][11,1]}

## IUCN90 
Coral_Current_IUCN90.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Current_IUCN90.mrl[k] <- Coral_PD.mrl[[k]][4,1]}
Coral_Past_IUCN90.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Past_IUCN90.mrl[k] <- Coral_PD.mrl[[k]][12,1]}

## IUCN100 
Coral_Current_IUCN100.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Current_IUCN100.mrl[k] <- Coral_PD.mrl[[k]][5,1]}
Coral_Past_IUCN100.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Past_IUCN100.mrl[k] <- Coral_PD.mrl[[k]][13,1]}

## IUCN490 
Coral_Current_IUCN490.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Current_IUCN490.mrl[k] <- Coral_PD.mrl[[k]][6,1]}
Coral_Past_IUCN490.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Past_IUCN490.mrl[k] <- Coral_PD.mrl[[k]][14,1]}

## IUCN500 
Coral_Current_IUCN500.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Current_IUCN500.mrl[k] <- Coral_PD.mrl[[k]][7,1]}
Coral_Past_IUCN500.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Past_IUCN500.mrl[k] <- Coral_PD.mrl[[k]][15,1]}

## Pess 
Coral_Current_Pess.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Current_Pess.mrl[k] <- Coral_PD.mrl[[k]][8,1]}
Coral_Past_Pess.mrl <- numeric(1000)
for (k in 1:1000) {Coral_Past_Pess.mrl[k] <- Coral_PD.mrl[[k]][16,1]}

Coral_results.mrl <- cbind(Coral_Past_Issac.mrl, Coral_Current_Issac.mrl, Coral_Past_IUCN40.mrl,Coral_Current_IUCN40.mrl, Coral_Past_IUCN50.mrl, Coral_Current_IUCN50.mrl,
                           Coral_Past_IUCN90.mrl, Coral_Current_IUCN90.mrl,Coral_Past_IUCN100.mrl, Coral_Current_IUCN100.mrl,Coral_Past_IUCN490.mrl, Coral_Current_IUCN490.mrl,
                           Coral_Past_IUCN500.mrl, Coral_Current_IUCN500.mrl,Coral_Past_Pess.mrl, Coral_Current_Pess.mrl)
Coral_results_stat.mrl <- stat.desc(Coral_results.mrl)


t1.1 <- t.test(Coral_Past_Issac,Coral_Past_Issac.mrl, paired = T)
t1.2 <- t.test(Coral_Current_Issac,Coral_Current_Issac.mrl, paired = T)

t1.3 <- t.test(Coral_Past_IUCN50,Coral_Past_IUCN50.mrl, paired = T) #
t1.4 <- t.test(Coral_Current_IUCN50,Coral_Current_IUCN50.mrl, paired = T)

t1.5 <- t.test(Coral_Past_IUCN100,Coral_Past_IUCN100.mrl, paired = T)
t1.6 <- t.test(Coral_Current_IUCN100,Coral_Current_IUCN100.mrl, paired = T)

t1.7 <- t.test(Coral_Past_IUCN500,Coral_Past_IUCN500.mrl, paired = T)
t1.8 <- t.test(Coral_Current_IUCN500,Coral_Current_IUCN500.mrl, paired = T)

t1.9 <- t.test(Coral_Past_Pess,Coral_Past_Pess.mrl, paired = T)
t1.10 <- t.test(Coral_Current_Pess,Coral_Current_Pess.mrl, paired = T)


g1.1 <- Hedgesg(Coral_Past_Issac,Coral_Past_Issac.mrl)
g1.2 <- Hedgesg(Coral_Current_Issac,Coral_Current_Issac.mrl)

g1.3 <- Hedgesg(Coral_Past_IUCN50,Coral_Past_IUCN50.mrl)
g1.4 <- Hedgesg(Coral_Current_IUCN50,Coral_Current_IUCN50.mrl)

g1.5 <- Hedgesg(Coral_Past_IUCN100,Coral_Past_IUCN100.mrl)
g1.6 <- Hedgesg(Coral_Current_IUCN100,Coral_Current_IUCN100.mrl)

g1.7 <- Hedgesg(Coral_Past_IUCN500,Coral_Past_IUCN500.mrl)
g1.8 <- Hedgesg(Coral_Current_IUCN500,Coral_Current_IUCN500.mrl)

g1.9 <- Hedgesg(Coral_Past_Pess,Coral_Past_Pess.mrl)
g1.10 <- Hedgesg(Coral_Current_Pess,Coral_Current_Pess.mrl)

g1.11 <- Hedgesg(Coral_Past_IUCN40,Coral_Past_IUCN40.mrl)
g1.12 <- Hedgesg(Coral_Current_IUCN40,Coral_Current_IUCN40.mrl)


t.all.pre <- cbind(t1.1,t1.2,t1.3,t1.4,t1.5,t1.6,t1.7,t1.8,t1.9,t1.10)
write.csv(t.all.pre, file = 'E:/you know files/IC/PD project/Data2/pre T test.csv')

g.all.pre <- cbind(g1.1,g1.2,g1.3,g1.4,g1.5,g1.6,g1.7,g1.8,g1.9,g1.10)
write.csv(g.all.pre, file = 'E:/you know files/IC/PD project/Data2/pre Hedgesg.csv')





##---------------------------------------------------------------------
## PD calculations --------------------

Coral_tree <- read.nexus(file = "E:/you know files/IC/PD project/Data2/Coral_tree_nosyn2.tre")

## PD calculation for 1000 trees
Coral_PD <- numeric(1000)
for (j in 1:1000) {
  Coral_Prob_input <- dataimpt(Coral_Prob_nosyn2)
  Coral_Prob_Num <- convertFormat(Coral_Prob_input)
  Coral_PD[j] <- ePD(tree = Coral_tree[[j]], probabilities.tips.present.matrix = Coral_Prob_Num)
}

## Results
## Issac
Coral_Current_Issac <- numeric(1000)
for (k in 1:1000) {Coral_Current_Issac[k] <- Coral_PD[[k]][1,1]}
Coral_Past_Issac <- numeric(1000)
for (k in 1:1000) {Coral_Past_Issac[k] <- Coral_PD[[k]][9,1]}

## IUCN40 
Coral_Current_IUCN40 <- numeric(1000)
for (k in 1:1000) {Coral_Current_IUCN40[k] <- Coral_PD[[k]][2,1]}
Coral_Past_IUCN40 <- numeric(1000)
for (k in 1:1000) {Coral_Past_IUCN40[k] <- Coral_PD[[k]][10,1]}

## IUCN50 
Coral_Current_IUCN50 <- numeric(1000)
for (k in 1:1000) {Coral_Current_IUCN50[k] <- Coral_PD[[k]][3,1]}
Coral_Past_IUCN50 <- numeric(1000)
for (k in 1:1000) {Coral_Past_IUCN50[k] <- Coral_PD[[k]][11,1]}

## IUCN90 
Coral_Current_IUCN90 <- numeric(1000)
for (k in 1:1000) {Coral_Current_IUCN90[k] <- Coral_PD[[k]][4,1]}
Coral_Past_IUCN90 <- numeric(1000)
for (k in 1:1000) {Coral_Past_IUCN90[k] <- Coral_PD[[k]][12,1]}

## IUCN100 
Coral_Current_IUCN100 <- numeric(1000)
for (k in 1:1000) {Coral_Current_IUCN100[k] <- Coral_PD[[k]][5,1]}
Coral_Past_IUCN100 <- numeric(1000)
for (k in 1:1000) {Coral_Past_IUCN100[k] <- Coral_PD[[k]][13,1]}

## IUCN490
Coral_Current_IUCN490 <- numeric(1000)
for (k in 1:1000) {Coral_Current_IUCN490[k] <- Coral_PD[[k]][6,1]}
Coral_Past_IUCN490 <- numeric(1000)
for (k in 1:1000) {Coral_Past_IUCN490[k] <- Coral_PD[[k]][14,1]}

## IUCN500 
Coral_Current_IUCN500 <- numeric(1000)
for (k in 1:1000) {Coral_Current_IUCN500[k] <- Coral_PD[[k]][7,1]}
Coral_Past_IUCN500 <- numeric(1000)
for (k in 1:1000) {Coral_Past_IUCN500[k] <- Coral_PD[[k]][15,1]}

## Pess 
Coral_Current_Pess <- numeric(1000)
for (k in 1:1000) {Coral_Current_Pess[k] <- Coral_PD[[k]][8,1]}
Coral_Past_Pess <- numeric(1000)
for (k in 1:1000) {Coral_Past_Pess[k] <- Coral_PD[[k]][16,1]}


Coral_results <- cbind(Coral_Past_Issac, Coral_Current_Issac, Coral_Past_IUCN40,Coral_Current_IUCN40, Coral_Past_IUCN50, Coral_Current_IUCN50,
                       Coral_Past_IUCN90, Coral_Current_IUCN90,Coral_Past_IUCN100, Coral_Current_IUCN100,Coral_Past_IUCN490, Coral_Current_IUCN490,
                       Coral_Past_IUCN500, Coral_Current_IUCN500,Coral_Past_Pess, Coral_Current_Pess)
Coral_results_stat <- stat.desc(Coral_results)
write.csv(Coral_results_stat, file = "E:/you know files/IC/PD project/Data2/Coral_results_stat.csv")  


##---------------------------------------------------------------------
## Randomization test------------------------------------

Coral_PD_random <- numeric(1000)
for (j in 1:1000) {
  Coral_Prob_input <- dataimpt(Coral_Prob_nosyn2)
  Coral_Prob_Num <- convertFormat(Coral_Prob_input)
  Coral_Prob_random <- Coral_Prob_Num
  rownames(Coral_Prob_random) <- sample(rownames(Coral_Prob_Num))
  Coral_PD_random[j] <- ePD(tree = Coral_tree[[j]], probabilities.tips.present.matrix = Coral_Prob_random)
}

## Issac
Random_Current_Issac <- numeric(1000)
for (k in 1:1000) {Random_Current_Issac[k] <- Coral_PD_random[[k]][1,1]}
Random_Past_Issac <- numeric(1000)
for (k in 1:1000) {Random_Past_Issac[k] <- Coral_PD_random[[k]][9,1]}

## IUCN40 
Random_Current_IUCN40 <- numeric(1000)
for (k in 1:1000) {Random_Current_IUCN40[k] <- Coral_PD_random[[k]][2,1]}
Random_Past_IUCN40 <- numeric(1000)
for (k in 1:1000) {Random_Past_IUCN40[k] <- Coral_PD_random[[k]][10,1]}

## IUCN50 
Random_Current_IUCN50 <- numeric(1000)
for (k in 1:1000) {Random_Current_IUCN50[k] <- Coral_PD_random[[k]][3,1]}
Random_Past_IUCN50 <- numeric(1000)
for (k in 1:1000) {Random_Past_IUCN50[k] <- Coral_PD_random[[k]][11,1]}

## IUCN90 
Random_Current_IUCN90 <- numeric(1000)
for (k in 1:1000) {Random_Current_IUCN90[k] <- Coral_PD_random[[k]][4,1]}
Random_Past_IUCN90 <- numeric(1000)
for (k in 1:1000) {Random_Past_IUCN90[k] <- Coral_PD_random[[k]][12,1]}

## IUCN100 
Random_Current_IUCN100 <- numeric(1000)
for (k in 1:1000) {Random_Current_IUCN100[k] <- Coral_PD_random[[k]][5,1]}
Random_Past_IUCN100 <- numeric(1000)
for (k in 1:1000) {Random_Past_IUCN100[k] <- Coral_PD_random[[k]][13,1]}

## IUCN490 
Random_Current_IUCN490 <- numeric(1000)
for (k in 1:1000) {Random_Current_IUCN490[k] <- Coral_PD_random[[k]][6,1]}
Random_Past_IUCN490 <- numeric(1000)
for (k in 1:1000) {Random_Past_IUCN490[k] <- Coral_PD_random[[k]][14,1]}

## IUCN500 
Random_Current_IUCN500 <- numeric(1000)
for (k in 1:1000) {Random_Current_IUCN500[k] <- Coral_PD_random[[k]][7,1]}
Random_Past_IUCN500 <- numeric(1000)
for (k in 1:1000) {Random_Past_IUCN500[k] <- Coral_PD_random[[k]][15,1]}

## Pess 
Random_Current_Pess <- numeric(1000)
for (k in 1:1000) {Random_Current_Pess[k] <- Coral_PD_random[[k]][8,1]}
Random_Past_Pess <- numeric(1000)
for (k in 1:1000) {Random_Past_Pess[k] <- Coral_PD_random[[k]][16,1]}

Random_results <- cbind(Random_Past_Issac, Random_Current_Issac, Random_Past_IUCN40,Random_Current_IUCN40, Random_Past_IUCN50, Random_Current_IUCN50,
                        Random_Past_IUCN90, Random_Current_IUCN90,Random_Past_IUCN100, Random_Current_IUCN100,Random_Past_IUCN490, Random_Current_IUCN490,
                        Random_Past_IUCN500, Random_Current_IUCN500,Random_Past_Pess, Random_Current_Pess)
Random_results_stat <- stat.desc(Random_results)
write.csv(Random_results_stat, file = "E:/you know files/IC/PD project/Data2/Random_results_stat.csv") 




## Sum of the TBL ---------------------
setwd("E:/you know files/IC/PD project/Data2")
library(dplyr)    ## dataframe manipulating
library(phytools) ## manipulating phylogenetic trees and data
library(pastecs) ## Package for Analysis of Space-Time Ecological Series (stat.desc Descriptive statistics on a data frame or time series)

Coral_Prob_binary <- select(Coral_Prob_nosyn2, 'Species', 'Current_Binary','Past_Binary')

Coral_Current_TBL <- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_binary_edited <- dataimptBi(Coral_Prob_binary)
  Current_threatened <- Coral_Prob_binary_edited%>%filter(Current_Binary == 1)
  Current_threatened.tip <- Current_threatened$Species
  
  # get the node number of the threatened tips
  Current_threatened.nodes<-sapply(Current_threatened.tip,function(x,y) which(y==x),y=Coral_tree_nosyn[[1]]$tip.label)
  
  # get the edge lengths for those nodes
  Current_edge.length<-setNames(Coral_tree_nosyn[[i]]$edge.length[sapply(Current_threatened.nodes,function(x,y) which(y==x),
                                                                         y=Coral_tree_nosyn[[i]]$edge[,2])],names(Current_threatened.nodes))
  Coral_Current_TBL[i] <- sum(Current_edge.length)
}

Coral_Past_TBL <- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_binary_edited <- dataimptBi(Coral_Prob_binary)
  Past_threatened <- Coral_Prob_binary_edited%>%filter(Past_Binary == 1)
  Past_threatened.tip <- Past_threatened$Species
  
  # get the node number of the threatened tips
  Past_threatened.nodes<-sapply(Past_threatened.tip,function(x,y) which(y==x),y=Coral_tree_nosyn[[1]]$tip.label)
  
  # get the edge lengths for those nodes
  Past_edge.length<-setNames(Coral_tree_nosyn[[i]]$edge.length[sapply(Past_threatened.nodes,function(x,y) which(y==x),
                                                                      y=Coral_tree_nosyn[[i]]$edge[,2])],names(Past_threatened.nodes))
  Coral_Past_TBL[i] <- sum(Past_edge.length)
}



## Randomization and randomized threatened edge length

Coral_Current_TBL.r <- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_binary_edited <- dataimptBi(Coral_Prob_binary)
  Current_threatened <- Coral_Prob_binary_edited%>%filter(Current_Binary == 1)
  
  Current_threatened.tip.r <- sample(Coral_Prob_binary_edited$Species, size = length(Current_threatened$Species))
  Current_threatened.nodes.r <- sapply(Current_threatened.tip.r,function(x,y) which(y==x),y=Coral_tree_nosyn[[i]]$tip.label)
  Current_edge.length.r <- setNames(Coral_tree_nosyn[[i]]$edge.length[sapply(Current_threatened.nodes.r, function(x,y) which(y==x),
                                                                             y=Coral_tree_nosyn[[i]]$edge[,2])],names(Current_threatened.nodes.r))
  Coral_Current_TBL.r[i] <- sum(Current_edge.length.r)
}


Coral_Past_TBL.r <- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_binary_edited <- dataimptBi(Coral_Prob_binary)
  Past_threatened <- Coral_Prob_binary_edited%>%filter(Past_Binary == 1)
  
  Past_threatened.tip.r <- sample(Coral_Prob_binary$Species, size = length(Past_threatened$Species))
  Past_threatened.nodes.r <- sapply(Past_threatened.tip.r,function(x,y) which(y==x),y=Coral_tree_nosyn[[i]]$tip.label)
  Past_edge.length.r <- setNames(Coral_tree_nosyn[[i]]$edge.length[sapply(Past_threatened.nodes.r,function(x,y) which(y==x),
                                                                          y=Coral_tree_nosyn[[i]]$edge[,2])],names(Past_threatened.nodes.r))
  Coral_Past_TBL.r[i] <- sum(Past_edge.length.r)
}

## Results
Coral_TBL <- cbind(Coral_Past_TBL, Coral_Past_TBL.r, Coral_Current_TBL,Coral_Current_TBL.r)
write.csv(Coral_TBL, "E:/you know files/IC/PD project/Data2/Coral_TBL.csv")

Coral_TBL.stat <- stat.desc(Coral_TBL)
write.csv(Coral_TBL.stat, "E:/you know files/IC/PD project/Data2/Coral_TBL_stat.csv")


g.Past.TBL <- Hedgesg(Coral_Past_TBL, Coral_Past_TBL.r)
g.Current.TBL <- Hedgesg(Coral_Current_TBL, Coral_Current_TBL.r)

g.TBL <- cbind(g.Past.TBL,g.Current.TBL)
write.csv(g.TBL, file = 'E:/you know files/IC/PD project/Data2/Hedgesg TBL.csv')

t13 <- t.test(Coral_Past_TBL, Coral_Past_TBL.r, paired = T)
t14 <- t.test(Coral_Current_TBL, Coral_Current_TBL.r, paired = T)

t.all.TBL <- cbind(t13,t14)

write.csv(t.all.TBL, file = 'E:/you know files/IC/PD project/Data2/t test TBL.csv')


## Binary PD and TBL---------------------------------
setwd("E:/you know files/IC/PD project/Data2")
library(dplyr)    ## dataframe manipulating
library(phytools) ## manipulating phylogenetic trees and data
library(pastecs) ## Package for Analysis of Space-Time Ecological Series (stat.desc Descriptive statistics on a data frame or time series)

Coral_Prob_binary <- select(Coral_Prob_nosyn2, 'Species', 'Current_Binary','Past_Binary')

## expected PD based on binary probability (pd1-pd2)
Coral_CurrentPD_binary <- numeric(1000)
for (l in 1:1000) {
  pd1 <- sum(Coral_tree_nosyn[[l]]$edge.length)
  Coral_Prob_binary_edited <- dataimptBi(Coral_Prob_binary)
  
  Current_threatened <- Coral_Prob_binary_edited%>%filter(Current_Binary == 1)
  Current_tree.new <- drop.tip(Coral_tree_nosyn[[l]], Coral_tree_nosyn[[l]]$tip.label[Coral_tree_nosyn[[l]]$tip.label %in% Current_threatened$Species])
  Current_pd2 <- sum(Current_tree.new$edge.length)
  Coral_CurrentPD_binary[l] <- pd1 - Current_pd2
}

Coral_PastPD_binary <- numeric(1000)
for (l in 1:1000) {
  pd1 <- sum(Coral_tree_nosyn[[l]]$edge.length)
  Coral_Prob_binary_edited <- dataimptBi(Coral_Prob_binary)
  
  Past_threatened <- filter(Coral_Prob_binary_edited, Past_Binary == 1)
  Past_tree.new <- drop.tip(Coral_tree_nosyn[[l]], Coral_tree_nosyn[[l]]$tip.label[Coral_tree_nosyn[[l]]$tip.label %in% Past_threatened$Species])
  Past_pd2 <- sum(Past_tree.new$edge.length)
  Coral_PastPD_binary[l] <- pd1 - Past_pd2
}


## randomization
Coral_CurrentPD_binary.r <- numeric(1000)
Coral_PastPD_binary.r <- numeric(1000)
for (i in 1:1000) {
  pd1 <- sum(Coral_tree_nosyn[[i]]$edge.length)
  Coral_Prob_binary_edited <- dataimptBi(Coral_Prob_binary)
  
  Current_threatened <- Coral_Prob_binary_edited%>%filter(Current_Binary == 1)
  Current_threatened.tip.r <- sample(Coral_Prob_binary_edited$Species, size = length(Current_threatened$Species))
  Current_tree.new.r <- drop.tip(Coral_tree_nosyn[[i]], Coral_tree_nosyn[[i]]$tip.label[Coral_tree_nosyn[[i]]$tip.label %in% Current_threatened.tip.r])
  Current_pd2.r <- sum(Current_tree.new.r$edge.length)
  Coral_CurrentPD_binary.r[i] <- pd1 - Current_pd2.r
  
  Past_threatened <- filter(Coral_Prob_binary_edited, Past_Binary == 1)
  Past_threatened.tip.r <- sample(Coral_Prob_binary_edited$Species, size = length(Past_threatened$Species))
  Past_tree.new.r <- drop.tip(Coral_tree_nosyn[[i]], Coral_tree_nosyn[[i]]$tip.label[Coral_tree_nosyn[[i]]$tip.label %in% Past_threatened.tip.r])
  Past_pd2.r <- sum(Past_tree.new.r$edge.length)
  Coral_PastPD_binary.r[i] <- pd1 - Past_pd2.r
}

#Results
Coral_PD_binary_both <- cbind(Coral_PastPD_binary,Coral_PastPD_binary.r,Coral_CurrentPD_binary,Coral_CurrentPD_binary.r)
write.csv(Coral_PD_binary_both, "E:/you know files/IC/PD project/Data2/pd1_pd2.csv")

Results_binary_both <- stat.desc(Coral_PD_binary_both)
write.csv(Results_binary_both, "E:/you know files/IC/PD project/Data2/Results_pd1_pd2.csv")

## Clustering------------------------------
install.packages('caper')
library(caper)
library(pastecs)

## Clustering analysis
D_Coral.2008 <- numeric(1000)
P.1.2008 <- numeric(1000)
P.0.2008 <- numeric(1000)
for (i in 1:1000) {
  Coral_Binary_edited <- dataimptBi(Coral_Prob_binary)
  rownames(Coral_Binary_edited) <- Coral_Binary_edited$Species
  Coral_Binary_edited$Current_Binary <- as.factor(Coral_Binary_edited$Current_Binary)
  D_test.2008 <- phylo.d(data = Coral_Binary_edited,phy = Coral_tree_nosyn[[i]],
                         names.col = Species,binvar = Current_Binary,permut = 1000)
  D_Coral.2008[i] <- D_test.2008$DEstimate
  P.1.2008[i] <- D_test.2008$Pval1
  P.0.2008[i] <- D_test.2008$Pval0
  print(i)
}

D_Coral.2008.stat <- stat.desc(D_Coral.2008)
P.1.2008.stat <- stat.desc(P.1.2008)
P.0.2008.stat <- stat.desc(P.0.2008)


D_Coral.1998 <- numeric(1000)
P.1.1998 <- numeric(1000)
P.0.1998 <- numeric(1000)
for (i in 1:1000) {
  Coral_Binary_edited <- dataimptBi(Coral_Prob_binary)
  rownames(Coral_Binary_edited) <- Coral_Binary_edited$Species
  Coral_Binary_edited$Past_Binary <- as.factor(Coral_Binary_edited$Past_Binary)
  D_test.1998 <- phylo.d(data = Coral_Binary_edited,phy = Coral_tree_nosyn[[i]],
                         names.col = Species,binvar = Past_Binary,permut = 1000)
  D_Coral.1998[i] <- D_test.1998$DEstimate
  P.1.1998[i] <- D_test.1998$Pval1
  P.0.1998[i] <- D_test.1998$Pval0
  print(i)
}

D_Coral.1998.stat <- stat.desc(D_Coral.1998)
P.1.1998.stat <- stat.desc(P.1.1998)
P.0.1998.stat <- stat.desc(P.0.1998)

D_Coral_stat <- cbind(D_Coral.2008.stat,P.1.2008.stat,P.0.2008.stat,
                      D_Coral.1998.stat,P.1.1998.stat,P.0.1998.stat)
write.csv(D_Coral_stat,file = "E:/you know files/IC/PD project/Data2/D_Coral_stat.csv")

D_Coral_all <- cbind(D_Coral.1998,P.1.1998,P.0.1998,D_Coral.2008,P.1.2008,P.0.2008)
write.csv(D_Coral_all,file = "E:/you know files/IC/PD project/Data2/D_Coral_stat.csv")



##---------------------------------------------------------------------
## percentage metrics------------------
Prob.ext.CR10 <-  0.5
Prob.ext.EN20 <-  0.2
Prob.ext.VU100 <-  0.1
Prob.ext.LC100 <- 0.0001
Prob.ext.NT100 <- 0.01

Prob.CR10 <- 1-Prob.ext.CR10
Prob.EN10 <- sqrt(1-Prob.ext.EN20)
Prob.VU10 <- (1-Prob.ext.VU100)^(1/10)
Prob.LC10 <- (1-Prob.ext.LC100)^(1/10)
Prob.NT10 <- (1-Prob.ext.NT100)^(1/10)

Prob.CR40 <- 1-Prob.CR10^4
Prob.CR50 <- 1-Prob.CR10^5
Prob.CR90 <- 1-Prob.CR10^9
Prob.CR100 <- 1-Prob.CR10^10
Prob.CR490 <- 1-Prob.CR10^49
Prob.CR500 <- 1-Prob.CR10^50

Prob.EN40 <- 1-Prob.EN10^4
Prob.EN50 <- 1-Prob.EN10^5
Prob.EN90 <- 1-Prob.EN10^9
Prob.EN100 <- 1-Prob.EN10^10
Prob.EN490 <- 1-Prob.EN10^49
Prob.EN500 <- 1-Prob.EN10^50

Prob.VU40 <- 1-Prob.VU10^4
Prob.VU50 <- 1-Prob.VU10^5
Prob.VU90 <- 1-Prob.VU10^9
Prob.VU100 <- 1-Prob.VU10^10
Prob.VU490 <- 1-Prob.VU10^49
Prob.VU500 <- 1-Prob.VU10^50

Prob.LC40 <- 1-Prob.LC10^4
Prob.LC50 <- 1-Prob.LC10^5
Prob.LC90 <- 1-Prob.LC10^9
Prob.LC100 <- 1-Prob.LC10^10
Prob.LC490 <- 1-Prob.LC10^49
Prob.LC500 <- 1-Prob.LC10^50

Prob.NT40 <- 1-Prob.NT10^4
Prob.NT50 <- 1-Prob.NT10^5
Prob.NT90 <- 1-Prob.NT10^9
Prob.NT100 <- 1-Prob.NT10^10
Prob.NT490 <- 1-Prob.NT10^49
Prob.NT500 <- 1-Prob.NT10^50

Prob.NT500.2 <- 1-(1-Prob.ext.NT100)^5

Coral_Prob_scaled <- matrix(nrow = 5, ncol = 6, byrow = FALSE, dimnames = list(c('LC', 'NT','VU','EN','CR'),
                                                                               c('IUCN40','IUCN50','IUCN90','IUCN100','IUCN490','IUCN500')))

Coral_Prob_scaled[5,] <- c(Prob.CR40,Prob.CR50,Prob.CR90,Prob.CR100,Prob.CR490,Prob.CR500)
Coral_Prob_scaled[4,] <- c(Prob.EN40,Prob.EN50,Prob.EN90,Prob.EN100,Prob.EN490,Prob.EN500)
Coral_Prob_scaled[3,] <- c(Prob.VU40,Prob.VU50,Prob.VU90,Prob.VU100,Prob.VU490,Prob.VU500)
Coral_Prob_scaled[2,] <- c(Prob.NT40,Prob.NT50,Prob.NT90,Prob.NT100,Prob.NT490,Prob.NT500)
Coral_Prob_scaled[1,] <- c(Prob.LC40,Prob.LC50,Prob.LC90,Prob.LC100,Prob.LC490,Prob.LC500)

write.csv(Coral_Prob_scaled, file = "E:/you know files/IC/PD project/Data/Coral_Prob_scaled.csv")

## % of PD change
per_change_Issac<- numeric(1000)
for (i in 1:1000) {
  per_change_Issac[i] <- (Coral_Current_Issac[i]-Coral_Past_Issac[i])/Coral_Past_Issac[i]*100
}

per_change_IUCN50<- numeric(1000)
for (i in 1:1000) {
  per_change_IUCN50[i] <- (Coral_Current_IUCN50[i]-Coral_Past_IUCN50[i])/Coral_Past_IUCN50[i]*100
}

per_change_IUCN100<- numeric(1000)
for (i in 1:1000) {
  per_change_IUCN100[i] <- (Coral_Current_IUCN100[i]-Coral_Past_IUCN100[i])/Coral_Past_IUCN100[i]*100
}

per_change_IUCN500<- numeric(1000)
for (i in 1:1000) {
  per_change_IUCN500[i] <- (Coral_Current_IUCN500[i]-Coral_Past_IUCN500[i])/Coral_Past_IUCN500[i]*100
}

per_change_Pess<- numeric(1000)
for (i in 1:1000) {
  per_change_Pess[i] <- (Coral_Current_Pess[i]-Coral_Past_Pess[i])/Coral_Past_Pess[i]*100
}

per_change_Binary<- numeric(1000)
for (i in 1:1000) {
  per_change_Binary[i] <- (Coral_CurrentPD_binary[i]-Coral_PastPD_binary[i])/Coral_PastPD_binary[i]*100
}

per_change <- cbind(per_change_Issac, per_change_IUCN50,per_change_IUCN100,per_change_IUCN500,per_change_Pess,per_change_Binary)
per_change.stat <- stat.desc(per_change)
write.csv(per_change.stat, file = "E:/you know files/IC/PD project/Data2/Coral_per_change.stat.csv")

# scaled % of change
per_change_IUCN50.s<- numeric(1000)
for (i in 1:1000) {
  per_change_IUCN50.s[i] <- (Coral_Current_IUCN40[i]-Coral_Past_IUCN50[i])/Coral_Past_IUCN50[i]*100
}

per_change_IUCN100.s<- numeric(1000)
for (i in 1:1000) {
  per_change_IUCN100.s[i] <- (Coral_Current_IUCN90[i]-Coral_Past_IUCN100[i])/Coral_Past_IUCN100[i]*100
}

per_change_IUCN500.s<- numeric(1000)
for (i in 1:1000) {
  per_change_IUCN500.s[i] <- (Coral_Current_IUCN490[i]-Coral_Past_IUCN500[i])/Coral_Past_IUCN500[i]*100
}

per_change.s <- cbind(per_change_IUCN50.s,per_change_IUCN100.s,per_change_IUCN500.s)
per_change.s.stat <- stat.desc(per_change.s)
write.csv(per_change.s.stat, file = "E:/you know files/IC/PD project/Data2/Coral_per_change_scaled.stat.csv")


## % of SR change
per_SR_Issac<- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_edited <- dataimpt(Coral_Prob)
  Coral_SR_Issac.1998 <- sum(Coral_Prob_edited$Past_Issac)
  Coral_SR_Issac.2008 <- sum(Coral_Prob_edited$Current_Issac)
  per_SR_Issac[i] <- (Coral_SR_Issac.2008-Coral_SR_Issac.1998)/Coral_SR_Issac.1998*100
}

per_SR_IUCN50<- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_edited <- dataimpt(Coral_Prob)
  Coral_SR_IUCN50.1998 <- sum(Coral_Prob_edited$Past_IUCN50)
  Coral_SR_IUCN50.2008 <- sum(Coral_Prob_edited$Current_IUCN50)
  per_SR_IUCN50[i] <- (Coral_SR_IUCN50.2008-Coral_SR_IUCN50.1998)/Coral_SR_IUCN50.1998*100
}

per_SR_IUCN100<- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_edited <- dataimpt(Coral_Prob)
  Coral_SR_IUCN100.1998 <- sum(Coral_Prob_edited$Past_IUCN100)
  Coral_SR_IUCN100.2008 <- sum(Coral_Prob_edited$Current_IUCN100)
  per_SR_IUCN100[i] <- (Coral_SR_IUCN100.2008-Coral_SR_IUCN100.1998)/Coral_SR_IUCN100.1998*100
}

per_SR_IUCN500<- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_edited <- dataimpt(Coral_Prob)
  Coral_SR_IUCN500.1998 <- sum(Coral_Prob_edited$Past_IUCN500)
  Coral_SR_IUCN500.2008 <- sum(Coral_Prob_edited$Current_IUCN500)
  per_SR_IUCN500[i] <- (Coral_SR_IUCN500.2008-Coral_SR_IUCN500.1998)/Coral_SR_IUCN500.1998*100
}

per_SR_Pess<- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_edited <- dataimpt(Coral_Prob)
  Coral_SR_Pess.1998 <- sum(Coral_Prob_edited$Past_Pess)
  Coral_SR_Pess.2008 <- sum(Coral_Prob_edited$Current_Pess)
  per_SR_Pess[i] <- (Coral_SR_Pess.2008-Coral_SR_Pess.1998)/Coral_SR_Pess.1998*100
}

per_SR_Binary<- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_edited <- dataimpt(Coral_Prob)
  Coral_SR_Binary.1998 <- sum(Coral_Prob_edited$Past_Binary)
  Coral_SR_Binary.2008 <- sum(Coral_Prob_edited$Current_Binary)
  per_SR_Binary[i] <- (Coral_SR_Binary.2008-Coral_SR_Binary.1998)/Coral_SR_Binary.1998*100
}

per_SR <- cbind(per_SR_Issac, per_SR_IUCN50,per_SR_IUCN100,per_SR_IUCN500,per_SR_Pess,per_SR_Binary)
per_SR.stat <- stat.desc(per_SR)
write.csv(per_SR.stat, file = "E:/you know files/IC/PD project/Data2/Coral_per_SR.stat.csv")

# scaled % of SR change
per_SR_IUCN50.s<- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_edited <- dataimpt(Coral_Prob)
  Coral_SR_IUCN50.1998 <- sum(Coral_Prob_edited$Past_IUCN50)
  Coral_SR_IUCN40.2008 <- sum(Coral_Prob_edited$Current_IUCN40)
  per_SR_IUCN50.s[i] <- (Coral_SR_IUCN40.2008-Coral_SR_IUCN50.1998)/Coral_SR_IUCN50.1998*100
}

per_SR_IUCN100.s<- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_edited <- dataimpt(Coral_Prob)
  Coral_SR_IUCN100.1998 <- sum(Coral_Prob_edited$Past_IUCN100)
  Coral_SR_IUCN90.2008 <- sum(Coral_Prob_edited$Current_IUCN90)
  per_SR_IUCN100.s[i] <- (Coral_SR_IUCN90.2008-Coral_SR_IUCN100.1998)/Coral_SR_IUCN100.1998*100
}

per_SR_IUCN500.s<- numeric(1000)
for (i in 1:1000) {
  Coral_Prob_edited <- dataimpt(Coral_Prob)
  Coral_SR_IUCN500.1998 <- sum(Coral_Prob_edited$Past_IUCN500)
  Coral_SR_IUCN490.2008 <- sum(Coral_Prob_edited$Current_IUCN490)
  per_SR_IUCN500.s[i] <- (Coral_SR_IUCN490.2008-Coral_SR_IUCN500.1998)/Coral_SR_IUCN500.1998*100
}

per_SR.s <- cbind(per_SR_IUCN50.s,per_SR_IUCN100.s,per_SR_IUCN500.s)
per_SR.s.stat <- stat.desc(per_SR.s)
write.csv(per_SR.s.stat, file = "E:/you know files/IC/PD project/Data2/Coral_per_SR_scaled.stat.csv")


##---------------------------------------------------------------------
## RLI and PDI---------------------------
##RLI script
wtEX<-5
wtPE<-5
wtEW<-5
wtPEW<-5
wtCR<-4
wtEN<-3
wtVU<-2
wtNT<-1
wtLC<-0

# Function to calculate the RLI from a set of species categories
get_rli <- function(species_categories) {
  # Get the number of species in each category, including the old categories
  nEX<-sum(species_categories=="EX" | species_categories=="E", na.rm = T)
  nPE<-sum(species_categories=="PE", na.rm = T)
  nEW<-sum(species_categories=="EW", na.rm = T)
  nPEW<-sum(species_categories=="PEW", na.rm = T)
  nCR<-sum(species_categories=="CR", na.rm = T)
  nEN<-sum(species_categories=="EN", na.rm = T)
  nVU<-sum(species_categories=="VU" | species_categories=="T", na.rm = T)
  nNT<-sum(species_categories=="NT"  | species_categories=="LR/nt", na.rm = T)
  nLC<-sum(species_categories=="LC" | species_categories=="LR/lc" | species_categories=="LR/cd", na.rm = T)
  
  # Multiply these by the respective weights and sum
  Total<-(wtEX*nEX)+(wtPE*nPE)+(wtEW*nEW)+(wtPEW*nPEW)+(wtCR*nCR)+(wtEN*nEN)+(wtVU*nVU)+(wtNT*nNT)+(wtLC*nLC)
  # Work out the worse case scenarios
  M<-((nEX+nPE+nEW+nPEW+nCR+nEN+nVU+nNT+nLC)*wtEX)
  # Calculate RLI 
  (RLI<-((M-Total)/M))
  return(RLI)
}

Coral_RL <- select(Coral_Prob_nosyn2, 'Species','Past_Cat','Current_Cat')
Coral_RLI.1998 <- numeric(1000)
Coral_RLI.2008 <- numeric(1000)
for (i in 1:1000) {
  Coral_RL_edited <- dataimptCat(Coral_RL)
  Coral_RLI.1998[i] <- get_rli(Coral_RL_edited$Past_Cat)
  Coral_RLI.2008[i] <- get_rli(Coral_RL_edited$Current_Cat)
}
Coral_RLI <- cbind(Coral_RLI.1998,Coral_RLI.2008)
Coral_RLI_stat <- cbind(stat.desc(Coral_RLI.1998),stat.desc(Coral_RLI.2008))
write.csv(Coral_RLI, file = "E:/you know files/IC/PD project/Data2/Coral_RLI.csv")
write.csv(Coral_RLI_stat, file = "E:/you know files/IC/PD project/Data2/Coral_RLI_stat.csv")

## PDI
tree <- list()
ePDloss <- numeric()
PDindex <- function(tree, ePDloss){
  PDI <- numeric(1000)
  for (i in 1:1000) {
    PD <- sum(tree[[i]]$edge.length)
    PDI[i] <- 1-ePDloss[i]/PD
  }
  return(PDI)
}
PDindex.m <- function(tree, ePDloss){
  PDI <- numeric(100)
  for (i in 1:100) {
    PD <- sum(tree[[i]]$edge.length)
    PDI[i] <- 1-ePDloss[i]/PD
  }
  return(PDI)
}

# Generating PDI for corals
Coral_PDI_Issac.1998 <- PDindex(Coral_tree_nosyn, Coral_Past_Issac)
Coral_PDI_Issac.2008 <- PDindex(Coral_tree_nosyn, Coral_Current_Issac)

Coral_PDI_IUCN40.1998 <- PDindex(Coral_tree_nosyn, Coral_Past_IUCN40)
Coral_PDI_IUCN40.2008 <- PDindex(Coral_tree_nosyn, Coral_Current_IUCN40)

Coral_PDI_IUCN50.1998 <- PDindex(Coral_tree_nosyn, Coral_Past_IUCN50)
Coral_PDI_IUCN50.2008 <- PDindex(Coral_tree_nosyn, Coral_Current_IUCN50)

Coral_PDI_IUCN90.1998 <- PDindex(Coral_tree_nosyn, Coral_Past_IUCN90)
Coral_PDI_IUCN90.2008 <- PDindex(Coral_tree_nosyn, Coral_Current_IUCN90)

Coral_PDI_IUCN100.1998 <- PDindex(Coral_tree_nosyn, Coral_Past_IUCN100)
Coral_PDI_IUCN100.2008 <- PDindex(Coral_tree_nosyn, Coral_Current_IUCN100)

Coral_PDI_IUCN490.1998 <- PDindex(Coral_tree_nosyn, Coral_Past_IUCN490)
Coral_PDI_IUCN490.2008 <- PDindex(Coral_tree_nosyn, Coral_Current_IUCN490)

Coral_PDI_IUCN500.1998 <- PDindex(Coral_tree_nosyn, Coral_Past_IUCN500)
Coral_PDI_IUCN500.2008 <- PDindex(Coral_tree_nosyn, Coral_Current_IUCN500)

Coral_PDI_Pess.1998 <- PDindex(Coral_tree_nosyn, Coral_Past_Pess)
Coral_PDI_Pess.2008 <- PDindex(Coral_tree_nosyn, Coral_Current_Pess)

Coral_PDI_Binary.1998 <- PDindex(Coral_tree_nosyn, Coral_PastPD_binary)
Coral_PDI_Binary.2008 <- PDindex(Coral_tree_nosyn, Coral_CurrentPD_binary)

Coral_PDI <- cbind(Coral_PDI_Issac.1998,Coral_PDI_Issac.2008,Coral_PDI_IUCN40.1998,Coral_PDI_IUCN40.2008,
                   Coral_PDI_IUCN50.1998,Coral_PDI_IUCN50.2008,Coral_PDI_IUCN90.1998,Coral_PDI_IUCN90.2008,
                   Coral_PDI_IUCN100.1998,Coral_PDI_IUCN100.2008,Coral_PDI_IUCN490.1998,Coral_PDI_IUCN490.2008,
                   Coral_PDI_IUCN500.1998,Coral_PDI_IUCN500.2008,Coral_PDI_Pess.1998,Coral_PDI_Pess.2008,
                   Coral_PDI_Binary.1998,Coral_PDI_Binary.2008)
Coral_PDI.stat <- stat.desc(Coral_PDI)
write.csv(Coral_PDI, file = "E:/you know files/IC/PD project/Data2/Coral_PDI.csv")
write.csv(Coral_PDI.stat, file = "E:/you know files/IC/PD project/Data2/Coral_PDI_stat.csv")


#Randomized PDI
Coral_rPDI_Issac.1998 <- PDindex(Coral_tree_nosyn, Random_Past_Issac)
Coral_rPDI_Issac.2008 <- PDindex(Coral_tree_nosyn, Random_Current_Issac)

Coral_rPDI_IUCN40.1998 <- PDindex(Coral_tree_nosyn, Random_Past_IUCN40)
Coral_rPDI_IUCN40.2008 <- PDindex(Coral_tree_nosyn, Random_Current_IUCN40)

Coral_rPDI_IUCN50.1998 <- PDindex(Coral_tree_nosyn, Random_Past_IUCN50)
Coral_rPDI_IUCN50.2008 <- PDindex(Coral_tree_nosyn, Random_Current_IUCN50)

Coral_rPDI_IUCN90.1998 <- PDindex(Coral_tree_nosyn, Random_Past_IUCN90)
Coral_rPDI_IUCN90.2008 <- PDindex(Coral_tree_nosyn, Random_Current_IUCN90)

Coral_rPDI_IUCN100.1998 <- PDindex(Coral_tree_nosyn, Random_Past_IUCN100)
Coral_rPDI_IUCN100.2008 <- PDindex(Coral_tree_nosyn, Random_Current_IUCN100)

Coral_rPDI_IUCN490.1998 <- PDindex(Coral_tree_nosyn, Random_Past_IUCN490)
Coral_rPDI_IUCN490.2008 <- PDindex(Coral_tree_nosyn, Random_Current_IUCN490)

Coral_rPDI_IUCN500.1998 <- PDindex(Coral_tree_nosyn, Random_Past_IUCN500)
Coral_rPDI_IUCN500.2008 <- PDindex(Coral_tree_nosyn, Random_Current_IUCN500)

Coral_rPDI_Pess.1998 <- PDindex(Coral_tree_nosyn, Random_Past_Pess)
Coral_rPDI_Pess.2008 <- PDindex(Coral_tree_nosyn, Random_Current_Pess)

Coral_rPDI_Binary.1998 <- PDindex(Coral_tree_nosyn, Coral_PastPD_binary.r)
Coral_rPDI_Binary.2008 <- PDindex(Coral_tree_nosyn, Coral_CurrentPD_binary.r)

Coral_rPDI <- cbind(Coral_rPDI_Issac.1998,Coral_rPDI_Issac.2008,Coral_rPDI_IUCN40.1998,Coral_rPDI_IUCN40.2008,
                    Coral_rPDI_IUCN50.1998,Coral_rPDI_IUCN50.2008,Coral_rPDI_IUCN90.1998,Coral_rPDI_IUCN90.2008,
                    Coral_rPDI_IUCN100.1998,Coral_rPDI_IUCN100.2008,Coral_rPDI_IUCN490.1998,Coral_rPDI_IUCN490.2008,
                    Coral_rPDI_IUCN500.1998,Coral_rPDI_IUCN500.2008,Coral_rPDI_Pess.1998,Coral_rPDI_Pess.2008,
                    Coral_rPDI_Binary.1998,Coral_rPDI_Binary.2008)
Coral_rPDI.stat <- stat.desc(Coral_rPDI)
write.csv(Coral_rPDI, file = "E:/you know files/IC/PD project/Data2/Coral_PDI_random.csv")
write.csv(Coral_rPDI.stat, file = "E:/you know files/IC/PD project/Data2/Coral_PDI_random_stat.csv")


##---------------------------------------------------------------------
## All of the codes for Mammals---------------------------------------------------
setwd("E:/you know files/IC/PD project/Data2/Mammal2")
library(devtools) # 开发小工具集合
library(phytools) ## manipulating phylogenetic trees and data
library(dplyr)    ## dataframe manipulating
library(ape)      ## phylogenetic analysis
library(mallorn)  ## PD calculation based on numeric prob of extinction
library(pastecs)  ## Package for Analysis of Space-Time Ecological Series (stat.desc Descriptive statistics on a data frame or time series)
library(ggplot2)

## preparing functions

convertProb.m <- function(Mam_RL){
  Mam_RL %>%
    mutate(Current_Issac = recode(Mam_RL$Current_Cat, 'LC' = 0.025, 'NT' = 0.05, 'VU' = 0.1, 'EN' = 0.2, 'CR' = 0.4,'EX' = 1, 'EW' = 1))%>%
    mutate(Current_IUCN38 = recode(Mam_RL$Current_Cat, 'LC' = 0.00004, 'NT' = 0.00381, 'VU' = 0.03925, 'EN' = 0.34556, 'CR' = 0.92821,'EX' = 1, 'EW' = 1)) %>%
    mutate(Current_IUCN50 = recode(Mam_RL$Current_Cat, 'LC' = 0.00005, 'NT' = 0.00501, 'VU' = 0.05132, 'EN' = 0.42757, 'CR' = 0.96875,'EX' = 1, 'EW' = 1)) %>%
    mutate(Current_IUCN88 = recode(Mam_RL$Current_Cat, 'LC' = 0.00009, 'NT' = 0.00881, 'VU' = 0.08855, 'EN' = 0.62538, 'CR' = 0.99776,'EX' = 1, 'EW' = 1))%>%
    mutate(Current_IUCN100 = recode(Mam_RL$Current_Cat, 'LC' = 0.00010, 'NT' = 0.01000, 'VU' = 0.10000, 'EN' = 0.67232, 'CR' = 0.99902,'EX' = 1, 'EW' = 1))%>%
    mutate(Current_IUCN488 = recode(Mam_RL$Current_Cat, 'LC' = 0.00049, 'NT' = 0.04786, 'VU' = 0.40200, 'EN' = 0.99568, 'CR' = 1,'EX' = 1, 'EW' = 1))%>%
    mutate(Current_IUCN500 = recode(Mam_RL$Current_Cat, 'LC' = 0.00050, 'NT' = 0.04901, 'VU' = 0.40951, 'EN' = 0.99622, 'CR' = 1,'EX' = 1, 'EW' = 1)) %>%
    mutate(Current_Pess = recode(Mam_RL$Current_Cat, 'LC' = 0.2, 'NT' = 0.4, 'VU' = 0.8, 'EN' = 0.9, 'CR' = 0.99,'EX' = 1, 'EW' = 1))%>%
    mutate(Current_Binary = recode(Mam_RL$Current_Cat, 'LC' = 0, 'NT' = 0, 'VU' = 1, 'EN' = 1, 'CR' = 1,'EX' = 1, 'EW' = 1))%>%
    
    mutate(Past_Issac = recode(Mam_RL$Past_Cat, 'LC' = 0.025, 'NT' = 0.05, 'VU' = 0.1, 'EN' = 0.2, 'CR' = 0.4,'EX' = 1, 'EW' = 1))%>%
    mutate(Past_IUCN38 = recode(Mam_RL$Past_Cat, 'LC' = 0.00004, 'NT' = 0.00381, 'VU' = 0.03925, 'EN' = 0.34556, 'CR' = 0.92821,'EX' = 1, 'EW' = 1)) %>%
    mutate(Past_IUCN50 = recode(Mam_RL$Past_Cat, 'LC' = 0.00005, 'NT' = 0.00501, 'VU' = 0.05132, 'EN' = 0.42757, 'CR' = 0.96875,'EX' = 1, 'EW' = 1)) %>%
    mutate(Past_IUCN88 = recode(Mam_RL$Past_Cat, 'LC' = 0.00009, 'NT' = 0.00881, 'VU' = 0.08855, 'EN' = 0.62538, 'CR' = 0.99776,'EX' = 1, 'EW' = 1))%>%
    mutate(Past_IUCN100 = recode(Mam_RL$Past_Cat, 'LC' = 0.00010, 'NT' = 0.01000, 'VU' = 0.10000, 'EN' = 0.67232, 'CR' = 0.99902,'EX' = 1, 'EW' = 1))%>%
    mutate(Past_IUCN488 = recode(Mam_RL$Past_Cat, 'LC' = 0.00049, 'NT' = 0.04786, 'VU' = 0.40200, 'EN' = 0.99568, 'CR' = 1,'EX' = 1, 'EW' = 1))%>%
    mutate(Past_IUCN500 = recode(Mam_RL$Past_Cat, 'LC' = 0.00050, 'NT' = 0.04901, 'VU' = 0.40951, 'EN' = 0.99622, 'CR' = 1,'EX' = 1, 'EW' = 1)) %>%
    mutate(Past_Pess = recode(Mam_RL$Past_Cat, 'LC' = 0.2, 'NT' = 0.4, 'VU' = 0.8, 'EN' = 0.9, 'CR' = 0.99,'EX' = 1, 'EW' = 1))%>%
    mutate(Past_Binary = recode(Mam_RL$Past_Cat, 'LC' = 0, 'NT' = 0, 'VU' = 1, 'EN' = 1, 'CR' = 1,'EX' = 1, 'EW' = 1))
}

dataimpt.m <- function(Coral_Prob_edited){
  Coral_Prob_na.index <- which(is.na(Coral_Prob_edited$Current_Issac)==T)     #取包含缺失值的行标
  Coral_Prob.clean <- Coral_Prob_edited[-Coral_Prob_na.index,]  
  Coral_Prob.na<-Coral_Prob_edited[Coral_Prob_na.index,]
  
  Coral_Prob.na$Current_Issac <-sample(Coral_Prob.clean$Current_Issac,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Current_IUCN38 <-sample(Coral_Prob.clean$Current_IUCN38,nrow(Coral_Prob.na),replace = T)
  Coral_Prob.na$Current_IUCN50 <-sample(Coral_Prob.clean$Current_IUCN50,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Current_IUCN88 <-sample(Coral_Prob.clean$Current_IUCN88,nrow(Coral_Prob.na),replace = T)
  Coral_Prob.na$Current_IUCN100 <-sample(Coral_Prob.clean$Current_IUCN100,nrow(Coral_Prob.na),replace = T)	
  Coral_Prob.na$Current_IUCN488 <-sample(Coral_Prob.clean$Current_IUCN488,nrow(Coral_Prob.na),replace = T)
  Coral_Prob.na$Current_IUCN500 <-sample(Coral_Prob.clean$Current_IUCN500,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Current_Pess <-sample(Coral_Prob.clean$Current_Pess,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Current_Binary <-sample(Coral_Prob.clean$Current_Binary,nrow(Coral_Prob.na),replace = T)		
  
  Coral_Prob.na$Past_Issac <-sample(Coral_Prob.clean$Past_Issac,nrow(Coral_Prob.na),replace = T)	
  Coral_Prob.na$Past_IUCN38<-sample(Coral_Prob.clean$Past_IUCN38,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_IUCN50<-sample(Coral_Prob.clean$Past_IUCN50,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_IUCN88<-sample(Coral_Prob.clean$Past_IUCN88,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_IUCN100 <-sample(Coral_Prob.clean$Past_IUCN100,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_IUCN488<-sample(Coral_Prob.clean$Past_IUCN488,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_IUCN500 <-sample(Coral_Prob.clean$Past_IUCN500,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_Pess <-sample(Coral_Prob.clean$Past_Pess,nrow(Coral_Prob.na),replace = T)		
  Coral_Prob.na$Past_Binary <-sample(Coral_Prob.clean$Past_Binary,nrow(Coral_Prob.na),replace = T)		
  
  rbind(Coral_Prob.clean, Coral_Prob.na)
}

convertFormat.m <- function(Prob_input){
  Prob_matrix <- select(Prob_input, -"Species",-"Current_Cat",-"Past_Cat",-"Current_Binary",-"Past_Binary")
  Prob_Num <- matrix(data = NA, nrow = dim(Prob_matrix)[1], ncol = dim(Prob_matrix)[2])
  for (i in 1:dim(Prob_matrix)[2]) {Prob_Num[,i] <- c(as.numeric(Prob_matrix[[i]]))}
  colnames(Prob_Num) <- c('Current_Issac','Current_IUCN38','Current_IUCN50','Current_IUCN88','Current_IUCN100','Current_IUCN488','Current_IUCN500','Current_Pess',
                          'Past_Issac','Past_IUCN38','Past_IUCN50','Past_IUCN88','Past_IUCN100','Past_IUCN488','Past_IUCN500','Past_Pess')
  rownames(Prob_Num) <- Prob_input$Species
  return(Prob_Num)
}

## PD calculation for 100 trees
Mam_tree <- read.nexus(file = "mammalphyloinclextinct.nex")
Mam_RL <- read.csv(file = "E:/you know files/IC/PD project/Data2/Mammal/Mam_RL.csv")
Mam_RL2 <- read.csv(file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_RL2.csv")
Mam_Prob2 <- convertProb.m(Mam_RL2)
write.csv(Mam_Prob2, file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_Prob2.csv")

EXEW <- filter(Mam_RL, Past_Cat == 'EX'|Past_Cat == 'EW')
tip.EXEW <- EXEW$Species
Mam_tree_noEX <- lapply(Mam_tree, drop.tip,tip=tip.EXEW)
class(Mam_tree_noEX)<-"multiPhylo"
writeNexus(Mam_tree_noEX, file="E:/you know files/IC/PD project/Data2/Mammal2/Mam_tree_noEX.tre")

Mam_PD_noEX <- numeric(100)
for (j in 1:100) {
  Mam_Prob_input <- dataimpt.m(Mam_Prob2)
  Mam_Prob_Num <- convertFormat.m(Mam_Prob_input)
  Mam_PD_noEX[j] <- ePD(tree = Mam_tree_noEX[[j]], probabilities.tips.present.matrix = Mam_Prob_Num)
  print(j)
}

## Results
## Issac
Mam_Current_Issac_noEX <- numeric(100)
for (k in 1:100) {Mam_Current_Issac_noEX[k] <- Mam_PD_noEX[[k]][1,1]}
Mam_Past_Issac_noEX <- numeric(100)
for (k in 1:100) {Mam_Past_Issac_noEX[k] <- Mam_PD_noEX[[k]][9,1]}

## IUCN38 
Mam_Current_IUCN38_noEX <- numeric(100)
for (k in 1:100) {Mam_Current_IUCN38_noEX[k] <- Mam_PD_noEX[[k]][2,1]}
Mam_Past_IUCN38_noEX <- numeric(100)
for (k in 1:100) {Mam_Past_IUCN38_noEX[k] <- Mam_PD_noEX[[k]][10,1]}

## IUCN50 
Mam_Current_IUCN50_noEX <- numeric(100)
for (k in 1:100) {Mam_Current_IUCN50_noEX[k] <- Mam_PD_noEX[[k]][3,1]}
Mam_Past_IUCN50_noEX <- numeric(100)
for (k in 1:100) {Mam_Past_IUCN50_noEX[k] <- Mam_PD_noEX[[k]][11,1]}

## IUCN88 
Mam_Current_IUCN88_noEX <- numeric(100)
for (k in 1:100) {Mam_Current_IUCN88_noEX[k] <- Mam_PD_noEX[[k]][4,1]}
Mam_Past_IUCN88_noEX <- numeric(100)
for (k in 1:100) {Mam_Past_IUCN88_noEX[k] <- Mam_PD_noEX[[k]][12,1]}

## IUCN100 
Mam_Current_IUCN100_noEX <- numeric(100)
for (k in 1:100) {Mam_Current_IUCN100_noEX[k] <- Mam_PD_noEX[[k]][5,1]}
Mam_Past_IUCN100_noEX <- numeric(100)
for (k in 1:100) {Mam_Past_IUCN100_noEX[k] <- Mam_PD_noEX[[k]][13,1]}

## IUCN488 
Mam_Current_IUCN488_noEX <- numeric(100)
for (k in 1:100) {Mam_Current_IUCN488_noEX[k] <- Mam_PD_noEX[[k]][6,1]}
Mam_Past_IUCN488_noEX <- numeric(100)
for (k in 1:100) {Mam_Past_IUCN488_noEX[k] <- Mam_PD_noEX[[k]][14,1]}

## IUCN500 
Mam_Current_IUCN500_noEX <- numeric(100)
for (k in 1:100) {Mam_Current_IUCN500_noEX[k] <- Mam_PD_noEX[[k]][7,1]}
Mam_Past_IUCN500_noEX <- numeric(100)
for (k in 1:100) {Mam_Past_IUCN500_noEX[k] <- Mam_PD_noEX[[k]][15,1]}

## Pess 
Mam_Current_Pess_noEX <- numeric(100)
for (k in 1:100) {Mam_Current_Pess_noEX[k] <- Mam_PD_noEX[[k]][8,1]}
Mam_Past_Pess_noEX <- numeric(100)
for (k in 1:100) {Mam_Past_Pess_noEX[k] <- Mam_PD_noEX[[k]][16,1]}

Mam_results_noEX <- cbind(Mam_Past_Issac_noEX, Mam_Current_Issac_noEX, Mam_Past_IUCN38_noEX,Mam_Current_IUCN38_noEX, Mam_Past_IUCN50_noEX, Mam_Current_IUCN50_noEX,
                          Mam_Past_IUCN88_noEX, Mam_Current_IUCN88_noEX,Mam_Past_IUCN100_noEX, Mam_Current_IUCN100_noEX,Mam_Past_IUCN488_noEX, Mam_Current_IUCN488_noEX,
                          Mam_Past_IUCN500_noEX, Mam_Current_IUCN500_noEX,Mam_Past_Pess_noEX, Mam_Current_Pess_noEX)
Mam_results_stat_noEX <- stat.desc(Mam_results_noEX)
write.csv(Mam_results_stat_noEX, file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_results_stat_noEX.csv")  

## randomised PD
Mam_PD_random_noEX <- numeric(100)
for (j in 1:100) {
  Mam_Prob_input <- dataimpt.m(Mam_Prob2)
  Mam_Prob_Num <- convertFormat.m(Mam_Prob_input)
  Mam_Prob_random <- Mam_Prob_Num
  rownames(Mam_Prob_random) <- sample(rownames(Mam_Prob_Num))
  Mam_PD_random_noEX[j] <- ePD(tree = Mam_tree_noEX[[j]], probabilities.tips.present.matrix = Mam_Prob_random)
  print(j)
}

## Result
## Issac
Random_Current_Issac.m_noEX <- numeric(100)
for (k in 1:100) {Random_Current_Issac.m_noEX[k] <- Mam_PD_random_noEX[[k]][1,1]}
Random_Past_Issac.m_noEX <- numeric(100)
for (k in 1:100) {Random_Past_Issac.m_noEX[k] <- Mam_PD_random_noEX[[k]][9,1]}

## IUCN38 
Random_Current_IUCN38.m_noEX <- numeric(100)
for (k in 1:100) {Random_Current_IUCN38.m_noEX[k] <- Mam_PD_random_noEX[[k]][2,1]}
Random_Past_IUCN38.m_noEX <- numeric(100)
for (k in 1:100) {Random_Past_IUCN38.m_noEX[k] <- Mam_PD_random_noEX[[k]][10,1]}

## IUCN50
Random_Current_IUCN50.m_noEX <- numeric(100)
for (k in 1:100) {Random_Current_IUCN50.m_noEX[k] <- Mam_PD_random_noEX[[k]][3,1]}
Random_Past_IUCN50.m_noEX <- numeric(100)
for (k in 1:100) {Random_Past_IUCN50.m_noEX[k] <- Mam_PD_random_noEX[[k]][11,1]}

## IUCN88
Random_Current_IUCN88.m_noEX <- numeric(100)
for (k in 1:100) {Random_Current_IUCN88.m_noEX[k] <- Mam_PD_random_noEX[[k]][4,1]}
Random_Past_IUCN88.m_noEX <- numeric(100)
for (k in 1:100) {Random_Past_IUCN88.m_noEX[k] <- Mam_PD_random_noEX[[k]][12,1]}

## IUCN100 
Random_Current_IUCN100.m_noEX <- numeric(100)
for (k in 1:100) {Random_Current_IUCN100.m_noEX[k] <- Mam_PD_random_noEX[[k]][5,1]}
Random_Past_IUCN100.m_noEX <- numeric(100)
for (k in 1:100) {Random_Past_IUCN100.m_noEX[k] <- Mam_PD_random_noEX[[k]][13,1]}

## IUCN488 
Random_Current_IUCN488.m_noEX <- numeric(100)
for (k in 1:100) {Random_Current_IUCN488.m_noEX[k] <- Mam_PD_random_noEX[[k]][6,1]}
Random_Past_IUCN488.m_noEX <- numeric(100)
for (k in 1:100) {Random_Past_IUCN488.m_noEX[k] <- Mam_PD_random_noEX[[k]][14,1]}

## IUCN500 
Random_Current_IUCN500.m_noEX <- numeric(100)
for (k in 1:100) {Random_Current_IUCN500.m_noEX[k] <- Mam_PD_random_noEX[[k]][7,1]}
Random_Past_IUCN500.m_noEX <- numeric(100)
for (k in 1:100) {Random_Past_IUCN500.m_noEX[k] <- Mam_PD_random_noEX[[k]][15,1]}

## Pess 
Random_Current_Pess.m_noEX <- numeric(100)
for (k in 1:100) {Random_Current_Pess.m_noEX[k] <- Mam_PD_random_noEX[[k]][8,1]}
Random_Past_Pess.m_noEX <- numeric(100)
for (k in 1:100) {Random_Past_Pess.m_noEX[k] <- Mam_PD_random_noEX[[k]][16,1]}

Random_results.m_noEX <- cbind(Random_Past_Issac.m_noEX, Random_Current_Issac.m_noEX, Random_Past_IUCN38.m_noEX,Random_Current_IUCN38.m_noEX, Random_Past_IUCN50.m_noEX, Random_Current_IUCN50.m_noEX,
                               Random_Past_IUCN88.m_noEX, Random_Current_IUCN88.m_noEX,Random_Past_IUCN100.m_noEX, Random_Current_IUCN100.m_noEX,Random_Past_IUCN488.m_noEX, Random_Current_IUCN488.m_noEX,
                               Random_Past_IUCN500.m_noEX, Random_Current_IUCN500.m_noEX,Random_Past_Pess.m_noEX, Random_Current_Pess.m_noEX)
Random_results_stat.m_noEX <- stat.desc(Random_results.m_noEX)
write.csv(Random_results_stat.m_noEX, file = "E:/you know files/IC/PD project/Data2/Mammal2/Random_results_stat_Mammal_noEX.csv") 



## expected PD based on binary probability (pd1-pd2)
Mam_Prob_binary_noEX <- select(Mam_Prob2, 'Species', 'Current_Binary','Past_Binary')
Mam_CurrentPD_binary_noEX <- numeric(100)
for (l in 1:100) {
  pd1 <- sum(Mam_tree_noEX[[l]]$edge.length)
  Mam_Prob_binary_edited <- dataimptBi(Mam_Prob_binary_noEX)
  
  Current_threatened <- Mam_Prob_binary_edited%>%filter(Current_Binary == 1)
  Current_tree.new <- drop.tip(Mam_tree_noEX[[l]], Mam_tree_noEX[[l]]$tip.label[Mam_tree_noEX[[l]]$tip.label %in% Current_threatened$Species])
  Current_pd2 <- sum(Current_tree.new$edge.length)
  Mam_CurrentPD_binary_noEX[l] <- pd1 - Current_pd2
}

Mam_PastPD_binary_noEX <- numeric(100)
for (l in 1:100) {
  pd1 <- sum(Mam_tree_noEX[[l]]$edge.length)
  Mam_Prob_binary_edited <- dataimptBi(Mam_Prob_binary_noEX)
  
  Past_threatened <- filter(Mam_Prob_binary_edited, Past_Binary == 1)
  Past_tree.new <- drop.tip(Mam_tree_noEX[[l]], Mam_tree_noEX[[l]]$tip.label[Mam_tree_noEX[[l]]$tip.label %in% Past_threatened$Species])
  Past_pd2 <- sum(Past_tree.new$edge.length)
  Mam_PastPD_binary_noEX[l] <- pd1 - Past_pd2
}


## binary randomization
Mam_CurrentPD_binary.r_noEX <- numeric(100)
Mam_PastPD_binary.r_noEX <- numeric(100)
for (i in 1:100) {
  pd1 <- sum(Mam_tree_noEX[[i]]$edge.length)
  Mam_Prob_binary_edited <- dataimptBi(Mam_Prob_binary_noEX)
  
  Current_threatened <- Mam_Prob_binary_edited%>%filter(Current_Binary == 1)
  Current_threatened.tip.r <- sample(Mam_Prob_binary_edited$Species, size = length(Current_threatened$Species))
  Current_tree.new.r <- drop.tip(Mam_tree_noEX[[i]], Mam_tree_noEX[[i]]$tip.label[Mam_tree_noEX[[i]]$tip.label %in% Current_threatened.tip.r])
  Current_pd2.r <- sum(Current_tree.new.r$edge.length)
  Mam_CurrentPD_binary.r_noEX[i] <- pd1 - Current_pd2.r
  
  Past_threatened <- filter(Mam_Prob_binary_edited, Past_Binary == 1)
  Past_threatened.tip.r <- sample(Mam_Prob_binary_edited$Species, size = length(Past_threatened$Species))
  Past_tree.new.r <- drop.tip(Mam_tree_noEX[[i]], Mam_tree_noEX[[i]]$tip.label[Mam_tree_noEX[[i]]$tip.label %in% Past_threatened.tip.r])
  Past_pd2.r <- sum(Past_tree.new.r$edge.length)
  Mam_PastPD_binary.r_noEX[i] <- pd1 - Past_pd2.r
}

#Results
Mam_PD_binary_both_noEX <- cbind(Mam_PastPD_binary_noEX,Mam_PastPD_binary.r_noEX,Mam_CurrentPD_binary_noEX,Mam_CurrentPD_binary.r_noEX)
write.csv(Mam_PD_binary_both_noEX, "E:/you know files/IC/PD project/Data2/Mammal2/pd1_pd2_noEX.csv")

Results_binary_both_noEX <- stat.desc(Mam_PD_binary_both_noEX)
write.csv(Results_binary_both_noEX, "E:/you know files/IC/PD project/Data2/Mammal2/Results_pd1_pd2_noEX.csv")


## RLI
Mam_RL_noEX <- select(Mam_Prob2, 'Species','Past_Cat','Current_Cat')
Mam_RLI.1996_noEX <- numeric(100)
Mam_RLI.2008_noEX <- numeric(100)
for (i in 1:100) {
  Mam_RL_edited <- dataimptCat(Mam_RL_noEX)
  Mam_RLI.1996_noEX[i] <- get_rli(Mam_RL_edited$Past_Cat)
  Mam_RLI.2008_noEX[i] <- get_rli(Mam_RL_edited$Current_Cat)
}
Mam_RLI_noEX <- cbind(Mam_RLI.1996_noEX,Mam_RLI.2008_noEX)
Mam_RLI_stat_noEX <- cbind(stat.desc(Mam_RLI.1996_noEX),stat.desc(Mam_RLI.2008_noEX))
write.csv(Mam_RLI_noEX, file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_RLI_noEX.csv")
write.csv(Mam_RLI_stat_noEX, file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_RLI_stat_noEX.csv")

## PDI
# Generating PDI for Mams
Mam_PDI_Issac.1996_noEX <- PDindex.m(Mam_tree_noEX , Mam_Past_Issac_noEX)
Mam_PDI_Issac.2008_noEX <- PDindex.m(Mam_tree_noEX , Mam_Current_Issac_noEX)

Mam_PDI_IUCN38.1996_noEX <- PDindex.m(Mam_tree_noEX , Mam_Past_IUCN38_noEX)
Mam_PDI_IUCN38.2008_noEX <- PDindex.m(Mam_tree_noEX , Mam_Current_IUCN38_noEX)

Mam_PDI_IUCN50.1996_noEX <- PDindex.m(Mam_tree_noEX , Mam_Past_IUCN50_noEX)
Mam_PDI_IUCN50.2008_noEX <- PDindex.m(Mam_tree_noEX , Mam_Current_IUCN50_noEX)

Mam_PDI_IUCN88.1996_noEX <- PDindex.m(Mam_tree_noEX, Mam_Past_IUCN88_noEX)
Mam_PDI_IUCN88.2008_noEX <- PDindex.m(Mam_tree_noEX, Mam_Current_IUCN88_noEX)

Mam_PDI_IUCN100.1996_noEX <- PDindex.m(Mam_tree_noEX, Mam_Past_IUCN100_noEX)
Mam_PDI_IUCN100.2008_noEX <- PDindex.m(Mam_tree_noEX, Mam_Current_IUCN100_noEX)

Mam_PDI_IUCN488.1996_noEX <- PDindex.m(Mam_tree_noEX, Mam_Past_IUCN488_noEX)
Mam_PDI_IUCN488.2008_noEX <- PDindex.m(Mam_tree_noEX, Mam_Current_IUCN488_noEX)

Mam_PDI_IUCN500.1996_noEX <- PDindex.m(Mam_tree_noEX, Mam_Past_IUCN500_noEX)
Mam_PDI_IUCN500.2008_noEX <- PDindex.m(Mam_tree_noEX, Mam_Current_IUCN500_noEX)

Mam_PDI_Pess.1996_noEX <- PDindex.m(Mam_tree_noEX, Mam_Past_Pess_noEX)
Mam_PDI_Pess.2008_noEX <- PDindex.m(Mam_tree_noEX, Mam_Current_Pess_noEX)

Mam_PDI_Binary.1996_noEX <- PDindex.m(Mam_tree_noEX, Mam_PastPD_binary_noEX)
Mam_PDI_Binary.2008_noEX <- PDindex.m(Mam_tree_noEX, Mam_CurrentPD_binary_noEX)

Mam_PDI_noEX <- cbind(Mam_PDI_Issac.1996_noEX,Mam_PDI_Issac.2008_noEX,Mam_PDI_IUCN38.1996_noEX,Mam_PDI_IUCN38.2008_noEX,
                      Mam_PDI_IUCN50.1996_noEX,Mam_PDI_IUCN50.2008_noEX,Mam_PDI_IUCN88.1996_noEX,Mam_PDI_IUCN88.2008_noEX,
                      Mam_PDI_IUCN100.1996_noEX,Mam_PDI_IUCN100.2008_noEX,Mam_PDI_IUCN488.1996_noEX,Mam_PDI_IUCN488.2008_noEX,
                      Mam_PDI_IUCN500.1996_noEX,Mam_PDI_IUCN500.2008_noEX,Mam_PDI_Pess.1996_noEX,Mam_PDI_Pess.2008_noEX,
                      Mam_PDI_Binary.1996_noEX,Mam_PDI_Binary.2008_noEX)
Mam_PDI.stat_noEX <- stat.desc(Mam_PDI_noEX)
write.csv(Mam_PDI_noEX, file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_PDI_noEX.csv")
write.csv(Mam_PDI.stat_noEX, file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_PDI_stat_noEX.csv")


#Randomized PDI
Mam_rPDI_Issac.1996_noEX <- PDindex.m(Mam_tree_noEX, Random_Past_Issac.m_noEX)
Mam_rPDI_Issac.2008_noEX <- PDindex.m(Mam_tree_noEX, Random_Current_Issac.m_noEX)

Mam_rPDI_IUCN38.1996_noEX <- PDindex.m(Mam_tree_noEX, Random_Past_IUCN38.m_noEX)
Mam_rPDI_IUCN38.2008_noEX <- PDindex.m(Mam_tree_noEX, Random_Current_IUCN38.m_noEX)

Mam_rPDI_IUCN50.1996_noEX <- PDindex.m(Mam_tree_noEX, Random_Past_IUCN50.m_noEX)
Mam_rPDI_IUCN50.2008_noEX <- PDindex.m(Mam_tree_noEX, Random_Current_IUCN50.m_noEX)

Mam_rPDI_IUCN88.1996_noEX <- PDindex.m(Mam_tree_noEX, Random_Past_IUCN88.m_noEX)
Mam_rPDI_IUCN88.2008_noEX <- PDindex.m(Mam_tree_noEX, Random_Current_IUCN88.m_noEX)

Mam_rPDI_IUCN100.1996_noEX <- PDindex.m(Mam_tree_noEX, Random_Past_IUCN100.m_noEX)
Mam_rPDI_IUCN100.2008_noEX <- PDindex.m(Mam_tree_noEX, Random_Current_IUCN100.m_noEX)

Mam_rPDI_IUCN488.1996_noEX <- PDindex.m(Mam_tree_noEX, Random_Past_IUCN488.m_noEX)
Mam_rPDI_IUCN488.2008_noEX <- PDindex.m(Mam_tree_noEX, Random_Current_IUCN488.m_noEX)

Mam_rPDI_IUCN500.1996_noEX <- PDindex.m(Mam_tree_noEX, Random_Past_IUCN500.m_noEX)
Mam_rPDI_IUCN500.2008_noEX <- PDindex.m(Mam_tree_noEX, Random_Current_IUCN500.m_noEX)

Mam_rPDI_Pess.1996_noEX <- PDindex.m(Mam_tree_noEX, Random_Past_Pess.m_noEX)
Mam_rPDI_Pess.2008_noEX <- PDindex.m(Mam_tree_noEX, Random_Current_Pess.m_noEX)

Mam_rPDI_Binary.1996_noEX <- PDindex.m(Mam_tree_noEX, Mam_PastPD_binary.r_noEX)
Mam_rPDI_Binary.2008_noEX <- PDindex.m(Mam_tree_noEX, Mam_CurrentPD_binary.r_noEX)

Mam_rPDI_noEX <- cbind(Mam_rPDI_Issac.1996_noEX,Mam_rPDI_Issac.2008_noEX,Mam_rPDI_IUCN38.1996_noEX,Mam_rPDI_IUCN38.2008_noEX,
                       Mam_rPDI_IUCN50.1996_noEX,Mam_rPDI_IUCN50.2008_noEX,Mam_rPDI_IUCN88.1996_noEX,Mam_rPDI_IUCN88.2008_noEX,
                       Mam_rPDI_IUCN100.1996_noEX,Mam_rPDI_IUCN100.2008_noEX,Mam_rPDI_IUCN488.1996_noEX,Mam_rPDI_IUCN488.2008_noEX,
                       Mam_rPDI_IUCN500.1996_noEX,Mam_rPDI_IUCN500.2008_noEX,Mam_rPDI_Pess.1996_noEX,Mam_rPDI_Pess.2008_noEX,
                       Mam_rPDI_Binary.1996_noEX,Mam_rPDI_Binary.2008_noEX)
Mam_rPDI.stat_noEX <- stat.desc(Mam_rPDI_noEX)
write.csv(Mam_rPDI_noEX, file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_PDI_random_noEX.csv")
write.csv(Mam_rPDI.stat_noEX, file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_PDI_random_stat_noEX.csv")



## % of change
Mam.per_change_Issac_noEX<- numeric(100)
for (i in 1:100) {
  Mam.per_change_Issac_noEX[i] <- (Mam_Current_Issac_noEX[i]-Mam_Past_Issac_noEX[i])/Mam_Past_Issac_noEX[i]*100
}

Mam.per_change_IUCN50_noEX<- numeric(100)
for (i in 1:100) {
  Mam.per_change_IUCN50_noEX[i] <- (Mam_Current_IUCN50_noEX[i]-Mam_Past_IUCN50_noEX[i])/Mam_Past_IUCN50_noEX[i]*100
}

Mam.per_change_IUCN100_noEX<- numeric(100)
for (i in 1:100) {
  Mam.per_change_IUCN100_noEX[i] <- (Mam_Current_IUCN100_noEX[i]-Mam_Past_IUCN100_noEX[i])/Mam_Past_IUCN100_noEX[i]*100
}

Mam.per_change_IUCN500_noEX<- numeric(100)
for (i in 1:100) {
  Mam.per_change_IUCN500_noEX[i] <- (Mam_Current_IUCN500_noEX[i]-Mam_Past_IUCN500_noEX[i])/Mam_Past_IUCN500_noEX[i]*100
}

Mam.per_change_Pess_noEX<- numeric(100)
for (i in 1:100) {
  Mam.per_change_Pess_noEX[i] <- (Mam_Current_Pess_noEX[i]-Mam_Past_Pess_noEX[i])/Mam_Past_Pess_noEX[i]*100
}

Mam.per_change_Binary_noEX<- numeric(100)
for (i in 1:100) {
  Mam.per_change_Binary_noEX[i] <- (Mam_CurrentPD_binary_noEX[i]-Mam_PastPD_binary_noEX[i])/Mam_PastPD_binary_noEX[i]*100
}

Mam.per_change_noEX <- cbind(Mam.per_change_Issac_noEX, Mam.per_change_IUCN50_noEX,Mam.per_change_IUCN100_noEX,Mam.per_change_IUCN500_noEX,Mam.per_change_Pess_noEX,Mam.per_change_Binary_noEX)
Mam.per_change.stat_noEX <- stat.desc(Mam.per_change_noEX)
write.csv(Mam.per_change.stat_noEX, file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_per_change.stat_noEX.csv")

# scaled % of change
Mam.per_change_IUCN50.s_noEX<- numeric(100)
for (i in 1:100) {
  Mam.per_change_IUCN50.s_noEX[i] <- (Mam_Current_IUCN38_noEX[i]-Mam_Past_IUCN50_noEX[i])/Mam_Past_IUCN50_noEX[i]*100
}

Mam.per_change_IUCN100.s_noEX<- numeric(100)
for (i in 1:100) {
  Mam.per_change_IUCN100.s_noEX[i] <- (Mam_Current_IUCN88_noEX[i]-Mam_Past_IUCN100_noEX[i])/Mam_Past_IUCN100_noEX[i]*100
}

Mam.per_change_IUCN500.s_noEX<- numeric(100)
for (i in 1:100) {
  Mam.per_change_IUCN500.s_noEX[i] <- (Mam_Current_IUCN488_noEX[i]-Mam_Past_IUCN500_noEX[i])/Mam_Past_IUCN500_noEX[i]*100
}

Mam.per_change.s_noEX <- cbind(Mam.per_change_IUCN50.s_noEX,Mam.per_change_IUCN100.s_noEX,Mam.per_change_IUCN500.s_noEX)
Mam.per_change.s.stat_noEX <- stat.desc(Mam.per_change.s_noEX)
write.csv(Mam.per_change.s.stat_noEX, file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_per_change_scaled.stat_noEX.csv")

## % of SR change
per_SR_Issac.m<- numeric(100)
for (i in 1:100) {
  Mam_Prob_edited <- dataimpt.m(Mam_Prob2)
  Mam_SR_Issac.1998 <- sum(Mam_Prob_edited$Past_Issac)
  Mam_SR_Issac.2008 <- sum(Mam_Prob_edited$Current_Issac)
  per_SR_Issac.m[i] <- (Mam_SR_Issac.2008-Mam_SR_Issac.1998)/Mam_SR_Issac.1998*100
}

per_SR_IUCN50.m<- numeric(100)
for (i in 1:100) {
  Mam_Prob_edited <- dataimpt.m(Mam_Prob2)
  Mam_SR_IUCN50.1998 <- sum(Mam_Prob_edited$Past_IUCN50)
  Mam_SR_IUCN50.2008 <- sum(Mam_Prob_edited$Current_IUCN50)
  per_SR_IUCN50.m[i] <- (Mam_SR_IUCN50.2008-Mam_SR_IUCN50.1998)/Mam_SR_IUCN50.1998*100
}

per_SR_IUCN100.m<- numeric(100)
for (i in 1:100) {
  Mam_Prob_edited <- dataimpt.m(Mam_Prob2)
  Mam_SR_IUCN100.1998 <- sum(Mam_Prob_edited$Past_IUCN100)
  Mam_SR_IUCN100.2008 <- sum(Mam_Prob_edited$Current_IUCN100)
  per_SR_IUCN100.m[i] <- (Mam_SR_IUCN100.2008-Mam_SR_IUCN100.1998)/Mam_SR_IUCN100.1998*100
}

per_SR_IUCN500.m<- numeric(100)
for (i in 1:100) {
  Mam_Prob_edited <- dataimpt.m(Mam_Prob2)
  Mam_SR_IUCN500.1998 <- sum(Mam_Prob_edited$Past_IUCN500)
  Mam_SR_IUCN500.2008 <- sum(Mam_Prob_edited$Current_IUCN500)
  per_SR_IUCN500.m[i] <- (Mam_SR_IUCN500.2008-Mam_SR_IUCN500.1998)/Mam_SR_IUCN500.1998*100
}

per_SR_Pess.m<- numeric(100)
for (i in 1:100) {
  Mam_Prob_edited <- dataimpt.m(Mam_Prob2)
  Mam_SR_Pess.1998 <- sum(Mam_Prob_edited$Past_Pess)
  Mam_SR_Pess.2008 <- sum(Mam_Prob_edited$Current_Pess)
  per_SR_Pess.m[i] <- (Mam_SR_Pess.2008-Mam_SR_Pess.1998)/Mam_SR_Pess.1998*100
}

per_SR_Binary.m<- numeric(100)
for (i in 1:100) {
  Mam_Prob_edited <- dataimpt.m(Mam_Prob2)
  Mam_SR_Binary.1998 <- sum(Mam_Prob_edited$Past_Binary)
  Mam_SR_Binary.2008 <- sum(Mam_Prob_edited$Current_Binary)
  per_SR_Binary.m[i] <- (Mam_SR_Binary.2008-Mam_SR_Binary.1998)/Mam_SR_Binary.1998*100
}

per_SR.m <- cbind(per_SR_Issac.m, per_SR_IUCN50.m,per_SR_IUCN100.m,per_SR_IUCN500.m,per_SR_Pess.m,per_SR_Binary.m)
per_SR.m.stat <- stat.desc(per_SR.m)
write.csv(per_SR.m.stat, file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_per_SR.stat.csv")

# scaled % of SR change
per_SR_IUCN50.m.s<- numeric(100)
for (i in 1:100) {
  Mam_Prob_edited <- dataimpt.m(Mam_Prob2)
  Mam_SR_IUCN50.1998 <- sum(Mam_Prob_edited$Past_IUCN50)
  Mam_SR_IUCN38.2008 <- sum(Mam_Prob_edited$Current_IUCN38)
  per_SR_IUCN50.m.s[i] <- (Mam_SR_IUCN38.2008-Mam_SR_IUCN50.1998)/Mam_SR_IUCN50.1998*100
}

per_SR_IUCN100.m.s<- numeric(100)
for (i in 1:100) {
  Mam_Prob_edited <- dataimpt.m(Mam_Prob2)
  Mam_SR_IUCN100.1998 <- sum(Mam_Prob_edited$Past_IUCN100)
  Mam_SR_IUCN88.2008 <- sum(Mam_Prob_edited$Current_IUCN88)
  per_SR_IUCN100.m.s[i] <- (Mam_SR_IUCN88.2008-Mam_SR_IUCN100.1998)/Mam_SR_IUCN100.1998*100
}

per_SR_IUCN500.m.s<- numeric(100)
for (i in 1:100) {
  Mam_Prob_edited <- dataimpt.m(Mam_Prob2)
  Mam_SR_IUCN500.1998 <- sum(Mam_Prob_edited$Past_IUCN500)
  Mam_SR_IUCN488.2008 <- sum(Mam_Prob_edited$Current_IUCN488)
  per_SR_IUCN500.m.s[i] <- (Mam_SR_IUCN488.2008-Mam_SR_IUCN500.1998)/Mam_SR_IUCN500.1998*100
}

per_SR.m.s <- cbind(per_SR_IUCN50.m.s,per_SR_IUCN100.m.s,per_SR_IUCN500.m.s)
per_SR.m.s.stat <- stat.desc(per_SR.m.s)
write.csv(per_SR.m.s.stat, file = "E:/you know files/IC/PD project/Data2/Mammal2/Mam_per_SR_scaled.stat.csv")



##---------------------------------------------------------------------
## All plot and test-----------------

library(pastecs) ## Package for Analysis of Space-Time Ecological Series (stat.desc Descriptive statistics on a data frame or time series)
install.packages('pwr')
library(pwr)

install.packages('ggpubr')
library(ggpubr)
## Hedge's g
g.Past.Issac <- Hedgesg(Coral_Past_Issac, Random_Past_Issac)
g.Past.IUCN50 <- Hedgesg(Coral_Past_IUCN50, Random_Past_IUCN50)
g.Past.IUCN100 <- Hedgesg(Coral_Past_IUCN100, Random_Past_IUCN100)
g.Past.IUCN500 <- Hedgesg(Coral_Past_IUCN500, Random_Past_IUCN500)
g.Past.Pess <- Hedgesg(Coral_Past_Pess, Random_Past_Pess)
g.Past.Binary <- Hedgesg(Coral_PastPD_binary, Coral_PastPD_binary.r)

g.Current.Issac <- Hedgesg(Coral_Current_Issac, Random_Current_Issac)
g.Current.IUCN50 <- Hedgesg(Coral_Current_IUCN50, Random_Current_IUCN50)
g.Current.IUCN100 <- Hedgesg(Coral_Current_IUCN100, Random_Current_IUCN100)
g.Current.IUCN500 <- Hedgesg(Coral_Current_IUCN500, Random_Current_IUCN500)
g.Current.Pess <- Hedgesg(Coral_Current_Pess, Random_Current_Pess)
g.Current.Binary <- Hedgesg(Coral_CurrentPD_binary, Coral_CurrentPD_binary.r)

g.ePDloss.random <- cbind(g.Past.Issac,g.Current.Issac,g.Past.IUCN50,g.Current.IUCN50,g.Past.IUCN100,g.Current.IUCN100,
                          g.Past.IUCN500,g.Current.IUCN500,g.Past.Pess,g.Current.Pess,g.Past.Binary,g.Current.Binary)
write.csv(g.ePDloss.random, file = 'E:/you know files/IC/PD project/Data2/Hedgesg random3.csv')



## Hedges'g for two time
g.Issac <- Hedgesg(Coral_Past_Issac, Coral_Current_Issac)
g.IUCN50 <- Hedgesg(Coral_Past_IUCN50, Coral_Current_IUCN50)
g.IUCN100 <- Hedgesg(Coral_Past_IUCN100, Coral_Current_IUCN100)
g.IUCN500 <- Hedgesg(Coral_Past_IUCN500, Coral_Current_IUCN500)
g.Pess <- Hedgesg(Coral_Past_Pess, Coral_Current_Pess)
g.Binary <- Hedgesg(Coral_PastPD_binary, Coral_CurrentPD_binary)

g.ePDloss.2time <- cbind(g.Issac,g.IUCN50,g.IUCN100,g.IUCN500,g.Pess,g.Binary)
write.csv(g.ePDloss.2time, file = 'E:/you know files/IC/PD project/Data2/Hedgesg 2time.csv')

## Plots
#Trends in expected PD loss through time

pdf('E:/you know files/IC/PD project/Data2/fig_trends_only3.pdf', width = 12, height = 8)  
par(mfrow = c(1,1), mar = c(8,6,6,6))
boxplot(Coral_Past_Issac , Coral_Current_Issac, 
        Coral_Past_IUCN50 ,Coral_Current_IUCN50, 
        Coral_Past_IUCN100 , Coral_Current_IUCN100, 
        Coral_Past_IUCN500 , Coral_Current_IUCN500, 
        Coral_Past_Pess ,Coral_Current_Pess, 
        Coral_PastPD_binary,Coral_CurrentPD_binary, 
        at = c(1,2, 4,5, 7,8, 10,11, 13,14, 16,17),
        col = c("lightslateblue","lightsalmon"),
        xaxt = 'n',
        ylim = c(0,6000)
)
abline(v = 3, col = 'gray72', lty = 2)
abline(v = 6, col = 'gray72', lty = 2)
abline(v = 9, col = 'gray72', lty = 2)
abline(v = 12, col = 'gray72', lty = 2)
abline(v = 15, col = 'gray72', lty = 2)

text(1.5,5800,"Issac", cex = 1.2)
text(4.5,5800,"IUCN50", cex = 1.2)
text(7.5,5800,"IUCN100", cex = 1.2)
text(10.5,5800,"IUCN500", cex = 1.2)
text(13.5,5800,"Pess.", cex = 1.2)
text(16.5,5800,"Binary", cex = 1.2)

mtext("1998", side=1, at=1, cex = 1.2, line = 0.5)
mtext("2008", side=1, at=2, cex = 1.2, line = 0.5)
mtext("1998", side=1, at=4, cex = 1.2, line = 0.5)
mtext("2008", side=1, at=5, cex = 1.2, line = 0.5)
mtext("1998", side=1, at=7, cex = 1.2, line = 0.5)
mtext("2008", side=1, at=8, cex = 1.2, line = 0.5)
mtext("1998", side=1, at=10, cex = 1.2, line = 0.5)
mtext("2008", side=1, at=11, cex = 1.2, line = 0.5)
mtext("1998", side=1, at=13, cex = 1.2, line = 0.5)
mtext("2008", side=1, at=14, cex = 1.2, line = 0.5)
mtext("1998", side=1, at=16, cex = 1.2, line = 0.5)
mtext("2008", side=1, at=17, cex = 1.2, line = 0.5)
mtext("Expected PD loss (Myr)", side = 2, line = 3, cex = 1.5)
mtext("Trends of Expected PD loss of reef-forming coral", side = 3, line = 2, font = 2, cex = 2)
legend(1,5000, legend = c('1998', '2008'), col = c("lightslateblue","lightsalmon"),
       pch = 15, cex = 1,bty = 'n')


dev.off()



## TBL
pdf('E:/you know files/IC/PD project/Data2/fig_TBL2.pdf', width = 7, height = 6)  
par(mar = c(4,8,4,5))
boxplot(Coral_Past_TBL, Coral_Past_TBL.r, Coral_Current_TBL,Coral_Current_TBL.r,
        at = c(1,2,4,5),
        names = c("1998", "1998", "2008", "2008"),
        col = c("lightgoldenrod1","olivedrab2"),
        las = 1,
        ylim = c(0, 1100))
legend(0.5,1000, legend = c('Threatened species', 'Random sample'), col = c("lightgoldenrod1","olivedrab2"),
       pch = 15, bty = 'n')
mtext("Sum of ternimal branch lengths (Myr)", side = 2, line = 3.5, font = 1, cex = 1.1 )
mtext("Expected TBL loss based on binary probability of extinction", side = 3, line = 2, font = 2, cex = 1.2 )
text(1.5,300,"Hedges' g  = -0.587")
text(4.5,1050,"Hedges' g  = -0.438")

dev.off()

## t-Test

t1 <- t.test(Coral_Past_Issac,Random_Past_Issac, paired = T)
t2 <- t.test(Coral_Current_Issac,Random_Current_Issac, paired = T)
pwr.t.test(d= g.Current.Issac,power = 0.8, sig.level = 0.05, type = 'two.sample' )

t3 <- t.test(Coral_Past_IUCN50,Random_Past_IUCN50, paired = T)
t4 <- t.test(Coral_Current_IUCN50,Random_Current_IUCN50, paired = T)

t5 <- t.test(Coral_Past_IUCN100,Random_Past_IUCN100, paired = T)
t6 <- t.test(Coral_Current_IUCN100,Random_Current_IUCN100, paired = T)

t7 <- t.test(Coral_Past_IUCN500,Random_Past_IUCN500, paired = T)
t8 <- t.test(Coral_Current_IUCN500,Random_Current_IUCN500, paired = T)

t9 <- t.test(Coral_Past_Pess,Random_Past_Pess, paired = T)
t10 <- t.test(Coral_Current_Pess,Random_Current_Pess, paired = T)

t11 <- t.test(Coral_PastPD_binary, Coral_PastPD_binary.r, paired = T)
t12 <- t.test(Coral_CurrentPD_binary, Coral_CurrentPD_binary.r, paired = T)


t.all.random <- cbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12)

write.csv(t.all.random, file = 'E:/you know files/IC/PD project/Data2/t test random.csv')


## t-test for 2 time points
t1.1 <- t.test(Coral_Past_Issac,Coral_Current_Issac, paired = T)

t3.1 <- t.test(Coral_Past_IUCN50,Coral_Current_IUCN50, paired = T)

t5.1 <- t.test(Coral_Past_IUCN100,Coral_Current_IUCN100, paired = T)

t7.1 <- t.test(Coral_Past_IUCN500,Coral_Current_IUCN500, paired = T)

t9.1 <- t.test(Coral_Past_Pess,Coral_Current_Pess, paired = T)

t11.1 <- t.test(Coral_PastPD_binary, Coral_CurrentPD_binary, paired = T)


t.all.2time <- cbind(t1.1,t3.1,t5.1,t7.1,t9.1,t11.1)

write.csv(t.all.2time, file = 'E:/you know files/IC/PD project/Data2/t test for 2 time points.csv')

## t-test for PDI and randomized PDI
t1.PDI.coral.Issac <- t.test(Coral_PDI_Issac.1998, Coral_rPDI_Issac.1998,paired = T)
t1.PDI.coral.IUCN50 <- t.test(Coral_PDI_IUCN50.1998, Coral_rPDI_IUCN50.1998,paired = T)
t1.PDI.coral.IUCN100 <- t.test(Coral_PDI_IUCN100.1998, Coral_rPDI_IUCN100.1998,paired = T)
t1.PDI.coral.IUCN500 <- t.test(Coral_PDI_IUCN500.1998, Coral_rPDI_IUCN50.1998,paired = T)
t1.PDI.coral.Pess <- t.test(Coral_PDI_Pess.1998, Coral_rPDI_Pess.1998,paired = T)
t1.PDI.coral.Binary <- t.test(Coral_PDI_Binary.1998, Coral_rPDI_Binary.1998,paired = T)

t2.PDI.coral.Issac <- t.test(Coral_PDI_Issac.2008, Coral_rPDI_Issac.2008,paired = T)
t2.PDI.coral.IUCN50 <- t.test(Coral_PDI_IUCN50.2008, Coral_rPDI_IUCN50.2008,paired = T)
t2.PDI.coral.IUCN100 <- t.test(Coral_PDI_IUCN100.2008, Coral_rPDI_IUCN100.2008,paired = T)
t2.PDI.coral.IUCN500 <- t.test(Coral_PDI_IUCN500.2008, Coral_rPDI_IUCN50.2008,paired = T)
t2.PDI.coral.Pess <- t.test(Coral_PDI_Pess.2008, Coral_rPDI_Pess.2008,paired = T)
t2.PDI.coral.Binary <- t.test(Coral_PDI_Binary.2008, Coral_rPDI_Binary.2008,paired = T)

t.PDI.coral <- cbind(t1.PDI.coral.Issac,t1.PDI.coral.IUCN50,t1.PDI.coral.IUCN100,t1.PDI.coral.IUCN500,
                     t1.PDI.coral.Pess,t1.PDI.coral.Binary,
                     t2.PDI.coral.Issac,t2.PDI.coral.IUCN50,t2.PDI.coral.IUCN100,t2.PDI.coral.IUCN500,
                     t2.PDI.coral.Pess,t2.PDI.coral.Binary)
write.csv(t.PDI.coral, file = 'E:/you know files/IC/PD project/Data2/t test PDI coral.csv')

g1.PDI.coral.Issac <- Hedgesg(Coral_PDI_Issac.1998, Coral_rPDI_Issac.1998)
g1.PDI.coral.IUCN50 <- Hedgesg(Coral_PDI_IUCN50.1998, Coral_rPDI_IUCN50.1998)
g1.PDI.coral.IUCN100 <- Hedgesg(Coral_PDI_IUCN100.1998, Coral_rPDI_IUCN100.1998)
g1.PDI.coral.IUCN500 <- Hedgesg(Coral_PDI_IUCN500.1998, Coral_rPDI_IUCN50.1998)
g1.PDI.coral.Pess <- Hedgesg(Coral_PDI_Pess.1998, Coral_rPDI_Pess.1998)
g1.PDI.coral.Binary <- Hedgesg(Coral_PDI_Binary.1998, Coral_rPDI_Binary.1998)

g2.PDI.coral.Issac <- Hedgesg(Coral_PDI_Issac.2008, Coral_rPDI_Issac.2008)
g2.PDI.coral.IUCN50 <- Hedgesg(Coral_PDI_IUCN50.2008, Coral_rPDI_IUCN50.2008)
g2.PDI.coral.IUCN100 <- Hedgesg(Coral_PDI_IUCN100.2008, Coral_rPDI_IUCN100.2008)
g2.PDI.coral.IUCN500 <- Hedgesg(Coral_PDI_IUCN500.2008, Coral_rPDI_IUCN50.2008)
g2.PDI.coral.Pess <- Hedgesg(Coral_PDI_Pess.2008, Coral_rPDI_Pess.2008)
g2.PDI.coral.Binary <- Hedgesg(Coral_PDI_Binary.2008, Coral_rPDI_Binary.2008)

g.PDI.coral <- cbind(g1.PDI.coral.Issac,g1.PDI.coral.IUCN50,g1.PDI.coral.IUCN100,g1.PDI.coral.IUCN500,
                     g1.PDI.coral.Pess,g1.PDI.coral.Binary,
                     g2.PDI.coral.Issac,g2.PDI.coral.IUCN50,g2.PDI.coral.IUCN100,g2.PDI.coral.IUCN500,
                     g2.PDI.coral.Pess,g2.PDI.coral.Binary)
write.csv(g.PDI.coral, file = 'E:/you know files/IC/PD project/Data2/Hedgesg PDI coral.csv')

t1.PDI.Mam.Issac <- t.test(Mam_PDI_Issac.1996_noEX, Mam_rPDI_Issac.1996_noEX,paired = T)
t1.PDI.Mam.IUCN50 <- t.test(Mam_PDI_IUCN50.1996_noEX, Mam_rPDI_IUCN50.1996_noEX,paired = T)
t1.PDI.Mam.IUCN100 <- t.test(Mam_PDI_IUCN100.1996_noEX, Mam_rPDI_IUCN100.1996_noEX,paired = T)
t1.PDI.Mam.IUCN500 <- t.test(Mam_PDI_IUCN500.1996_noEX, Mam_rPDI_IUCN50.1996_noEX,paired = T)
t1.PDI.Mam.Pess <- t.test(Mam_PDI_Pess.1996_noEX, Mam_rPDI_Pess.1996_noEX,paired = T)
t1.PDI.Mam.Binary <- t.test(Mam_PDI_Binary.1996_noEX, Mam_rPDI_Binary.1996_noEX,paired = T)

t2.PDI.Mam.Issac <- t.test(Mam_PDI_Issac.2008_noEX, Mam_rPDI_Issac.2008_noEX,paired = T)
t2.PDI.Mam.IUCN50 <- t.test(Mam_PDI_IUCN50.2008_noEX, Mam_rPDI_IUCN50.2008_noEX,paired = T)
t2.PDI.Mam.IUCN100 <- t.test(Mam_PDI_IUCN100.2008_noEX, Mam_rPDI_IUCN100.2008_noEX,paired = T)
t2.PDI.Mam.IUCN500 <- t.test(Mam_PDI_IUCN500.2008_noEX, Mam_rPDI_IUCN50.2008_noEX,paired = T)
t2.PDI.Mam.Pess <- t.test(Mam_PDI_Pess.2008_noEX, Mam_rPDI_Pess.2008_noEX,paired = T)
t2.PDI.Mam.Binary <- t.test(Mam_PDI_Binary.2008_noEX, Mam_rPDI_Binary.2008_noEX,paired = T)

t.PDI.Mam <- cbind(t1.PDI.Mam.Issac,t1.PDI.Mam.IUCN50,t1.PDI.Mam.IUCN100,t1.PDI.Mam.IUCN500,
                   t1.PDI.Mam.Pess,t1.PDI.Mam.Binary,
                   t2.PDI.Mam.Issac,t2.PDI.Mam.IUCN50,t2.PDI.Mam.IUCN100,t2.PDI.Mam.IUCN500,
                   t2.PDI.Mam.Pess,t2.PDI.Mam.Binary)
write.csv(t.PDI.Mam, file = 'E:/you know files/IC/PD project/Data2/t test PDI Mam.csv')

g1.PDI.Mam.Issac <- Hedgesg(Mam_PDI_Issac.1996_noEX, Mam_rPDI_Issac.1996_noEX)
g1.PDI.Mam.IUCN50 <- Hedgesg(Mam_PDI_IUCN50.1996_noEX, Mam_rPDI_IUCN50.1996_noEX)
g1.PDI.Mam.IUCN100 <- Hedgesg(Mam_PDI_IUCN100.1996_noEX, Mam_rPDI_IUCN100.1996_noEX)
g1.PDI.Mam.IUCN500 <- Hedgesg(Mam_PDI_IUCN500.1996_noEX, Mam_rPDI_IUCN50.1996_noEX)
g1.PDI.Mam.Pess <- Hedgesg(Mam_PDI_Pess.1996_noEX, Mam_rPDI_Pess.1996_noEX)
g1.PDI.Mam.Binary <- Hedgesg(Mam_PDI_Binary.1996_noEX, Mam_rPDI_Binary.1996_noEX)

g2.PDI.Mam.Issac <- Hedgesg(Mam_PDI_Issac.2008_noEX, Mam_rPDI_Issac.2008_noEX)
g2.PDI.Mam.IUCN50 <- Hedgesg(Mam_PDI_IUCN50.2008_noEX, Mam_rPDI_IUCN50.2008_noEX)
g2.PDI.Mam.IUCN100 <- Hedgesg(Mam_PDI_IUCN100.2008_noEX, Mam_rPDI_IUCN100.2008_noEX)
g2.PDI.Mam.IUCN500 <- Hedgesg(Mam_PDI_IUCN500.2008_noEX, Mam_rPDI_IUCN50.2008_noEX)
g2.PDI.Mam.Pess <- Hedgesg(Mam_PDI_Pess.2008_noEX, Mam_rPDI_Pess.2008_noEX)
g2.PDI.Mam.Binary <- Hedgesg(Mam_PDI_Binary.2008_noEX, Mam_rPDI_Binary.2008_noEX)

g.PDI.Mam <- cbind(g1.PDI.Mam.Issac,g1.PDI.Mam.IUCN50,g1.PDI.Mam.IUCN100,g1.PDI.Mam.IUCN500,
                   g1.PDI.Mam.Pess,g1.PDI.Mam.Binary,
                   g2.PDI.Mam.Issac,g2.PDI.Mam.IUCN50,g2.PDI.Mam.IUCN100,g2.PDI.Mam.IUCN500,
                   g2.PDI.Mam.Pess,g2.PDI.Mam.Binary)
write.csv(g.PDI.Mam, file = 'E:/you know files/IC/PD project/Data2/Hedgesg PDI Mam.csv')


library(ggplot2)
## PDI and RLI
year <- c(1998,2008)
Coral_RLI.mean <- c(0.985061945,0.808951501)
Coral_PDI_Issac.mean <- c(0.70555725,0.607724275)
Coral_PDI_IUCN50.mean <- c(0.836927752,0.682123385)
Coral_PDI_IUCN100.mean <- c(0.814827841,0.611934311)
Coral_PDI_IUCN500.mean <- c(0.779143818,0.433495292)
Coral_PDI_Pess.mean <- c(0.37159903,0.220086305)
Coral_PDI_Binary.mean <- c(0.996372637,0.889859182)

## Plots for RLI and randomized PDI -

year <- c(1998,2008)
Coral_RLI.mean <- c(0.985061945,0.808951501)
Coral_PDI_Issac.mean <- c(0.70555725,0.607724275)
Coral_PDI_IUCN50.mean <- c(0.836927752,0.682123385)
Coral_PDI_IUCN100.mean <- c(0.814827841,0.611934311)
Coral_PDI_IUCN500.mean <- c(0.779143818,0.433495292)
Coral_PDI_Pess.mean <- c(0.37159903,0.220086305)
Coral_PDI_Binary.mean <- c(0.996372637,0.889859182)

Coral_PDI_Issac.r.mean <- c(0.707629611,0.605151442)
Coral_PDI_IUCN50.r.mean <- c(0.839041879,0.678208967)
Coral_PDI_IUCN100.r.mean <- c(0.815823951,0.607571504)
Coral_PDI_IUCN500.r.mean <- c(0.772759313,0.422961585)
Coral_PDI_Pess.r.mean <- c(0.373862623,0.21065784)
Coral_PDI_Binary.r.mean <- c(0.994496854,0.888206866)

year.m <- c(1996,2008)
Mam_RLI.mean_noEX <- c(0.85791358,0.851319273)
Mam_PDI_Issac.mean_noEX <- c(0.79235303,0.786337489)
Mam_PDI_IUCN50.mean_noEX <- c(0.777354821,0.764179735)
Mam_PDI_IUCN100.mean_noEX <- c(0.731665922,0.719670125)
Mam_PDI_IUCN500.mean_noEX <- c(0.644643646,0.633478868)
Mam_PDI_Pess.mean_noEX <- c(0.425411427,0.419127489)
Mam_PDI_Binary.mean_noEX <- c(0.876290801,0.869055307)

Mam_PDI_Issac.r.mean_noEX <- c(0.787532995,0.781497516)
Mam_PDI_IUCN50.r.mean_noEX <- c( 0.758157871,0.744705187)
Mam_PDI_IUCN100.r.mean_noEX <- c(0.706811504,0.69390833)
Mam_PDI_IUCN500.r.mean_noEX <- c(0.612278323,0.600118796)
Mam_PDI_Pess.r.mean_noEX <- c(0.406903597,0.400846132)
Mam_PDI_Binary.r.mean_noEX <- c(0.872237529,0.86643087)


dd <- data.frame(year,Coral_RLI.mean,Coral_PDI_Issac.mean,Coral_PDI_IUCN50.mean,Coral_PDI_IUCN100.mean,
                 Coral_PDI_IUCN500.mean,Coral_PDI_Pess.mean,Coral_PDI_Binary.mean)
dd.c <- data.frame(year,Coral_RLI.mean,Coral_PDI_Issac.mean,Coral_PDI_IUCN50.mean,Coral_PDI_IUCN100.mean,
                   Coral_PDI_IUCN500.mean,Coral_PDI_Pess.mean,Coral_PDI_Binary.mean,
                   Coral_PDI_Issac.r.mean,Coral_PDI_IUCN50.r.mean,Coral_PDI_IUCN100.r.mean,
                   Coral_PDI_IUCN500.r.mean,Coral_PDI_Pess.r.mean,Coral_PDI_Binary.r.mean)
dd.m_noEX <- data.frame(year.m,Mam_RLI.mean_noEX,Mam_PDI_Issac.mean_noEX,Mam_PDI_IUCN50.mean_noEX,Mam_PDI_IUCN100.mean_noEX,
                        Mam_PDI_IUCN500.mean_noEX,Mam_PDI_Pess.mean_noEX,Mam_PDI_Binary.mean_noEX)
dd.m.c_noEX <- data.frame(year.m,Mam_RLI.mean_noEX,Mam_PDI_Issac.mean_noEX,Mam_PDI_IUCN50.mean_noEX,Mam_PDI_IUCN100.mean_noEX,
                          Mam_PDI_IUCN500.mean_noEX,Mam_PDI_Pess.mean_noEX,Mam_PDI_Binary.mean_noEX,
                          Mam_PDI_Issac.r.mean_noEX,Mam_PDI_IUCN50.r.mean_noEX,Mam_PDI_IUCN100.r.mean_noEX,
                          Mam_PDI_IUCN500.r.mean_noEX,Mam_PDI_Pess.r.mean_noEX,Mam_PDI_Binary.r.mean_noEX)


## combine PDI and randomized PDI of corals 

Plot_PDI_Coral.r <-ggplot()+ 
  geom_line(data = dd.c,aes(x = year,y = Coral_RLI.mean,colour = "RLI"),size=0.5)+ 
  geom_point(data = dd.c,aes(x = year,y = Coral_RLI.mean,colour = "RLI"),size=1)+
  geom_line(data = dd.c,aes(x = year,y = Coral_PDI_Issac.mean,colour = 'PDI_Issac',linetype="RLI or PDI"),size=0.5)+
  geom_point(data = dd.c,aes(x = year,y = Coral_PDI_Issac.mean,colour = 'PDI_Issac'),size=1)+
  geom_line(data = dd.c,aes(x = year,y = Coral_PDI_IUCN50.mean,colour = 'PDI_IUCN50',linetype="RLI or PDI"),size=0.5)+
  geom_point(data = dd.c,aes(x = year,y = Coral_PDI_IUCN50.mean,colour = 'PDI_IUCN50'),size=1)+ 
  geom_line(data = dd.c,aes(x = year,y = Coral_PDI_IUCN100.mean,colour = 'PDI_IUCN100',linetype="RLI or PDI"),size=0.5)+
  geom_point(data = dd.c,aes(x = year,y = Coral_PDI_IUCN100.mean,colour = 'PDI_IUCN100'),size=1)+
  geom_line(data = dd.c,aes(x = year,y = Coral_PDI_IUCN500.mean,colour = 'PDI_IUCN500',linetype="RLI or PDI"),size=0.5)+
  geom_point(data = dd.c,aes(x = year,y = Coral_PDI_IUCN500.mean,colour = 'PDI_IUCN500'),size=1)+
  geom_line(data = dd.c,aes(x = year,y = Coral_PDI_Pess.mean,colour = 'PDI_Pessimistic',linetype="RLI or PDI"),size=0.5)+
  geom_point(data = dd.c,aes(x = year,y = Coral_PDI_Pess.mean,colour = 'PDI_Pessimistic'),size=1)+
  geom_line(data = dd.c,aes(x = year,y = Coral_PDI_Binary.mean,colour = 'PDI_Binary',linetype="RLI or PDI"),size=0.5)+
  geom_point(data = dd.c,aes(x = year,y = Coral_PDI_Binary.mean,colour = 'PDI_Binary'),size=1)+
  
  geom_line(data = dd.c,aes(x = year,y = Coral_PDI_Issac.r.mean,colour = 'PDI_Issac.r',linetype="randomized PDI"),size=0.5)+
  geom_point(data = dd.c,aes(x = year,y = Coral_PDI_Issac.r.mean,colour = 'PDI_Issac.r'),size=1)+
  geom_line(data = dd.c,aes(x = year,y = Coral_PDI_IUCN50.r.mean,colour = 'PDI_IUCN50.r',linetype="randomized PDI"),size=0.5)+
  geom_point(data = dd.c,aes(x = year,y = Coral_PDI_IUCN50.r.mean,colour = 'PDI_IUCN50.r'),size=1)+ 
  geom_line(data = dd.c,aes(x = year,y = Coral_PDI_IUCN100.r.mean,colour = 'PDI_IUCN100.r',linetype="randomized PDI"),size=0.5)+
  geom_point(data = dd.c,aes(x = year,y = Coral_PDI_IUCN100.r.mean,colour = 'PDI_IUCN100.r'),size=1)+
  geom_line(data = dd.c,aes(x = year,y = Coral_PDI_IUCN500.r.mean,colour = 'PDI_IUCN500.r',linetype="randomized PDI"),size=0.5)+
  geom_point(data = dd.c,aes(x = year,y = Coral_PDI_IUCN500.r.mean,colour = 'PDI_IUCN500.r'),size=1)+
  geom_line(data = dd.c,aes(x = year,y = Coral_PDI_Pess.r.mean,colour = 'PDI_Pessimistic.r',linetype="randomized PDI"),size=0.5)+
  geom_point(data = dd.c,aes(x = year,y = Coral_PDI_Pess.r.mean,colour = 'PDI_Pessimistic.r'),size=1)+
  geom_line(data = dd.c,aes(x = year,y = Coral_PDI_Binary.r.mean,colour = 'PDI_Binary.r',linetype="randomized PDI"),size=0.5)+
  geom_point(data = dd.c,aes(x = year,y = Coral_PDI_Binary.r.mean,colour = 'PDI_Binary.r'),size=1)+
  
  ylim(0.2,1.1)+scale_x_continuous(breaks=seq(1996, 2008, 4))+
  scale_colour_manual('', values = c('RLI' = "Red",'PDI_Issac' = 'pink','PDI_IUCN50' = 'steelblue1','PDI_IUCN100' = 'steelblue3',
                                     'PDI_IUCN500' = 'steelblue4','PDI_Pessimistic' = 'plum3','PDI_Binary' = 'plum4','PDI_Issac.r' = 'pink',
                                     'PDI_IUCN50.r' = 'steelblue1','PDI_IUCN100.r' = 'steelblue3',
                                     'PDI_IUCN500.r' = 'steelblue4','PDI_Pessimistic.r' = 'plum3','PDI_Binary.r' = 'plum4'),
                      breaks = c('RLI','PDI_Issac','PDI_IUCN50','PDI_IUCN100','PDI_IUCN500','PDI_Pessimistic','PDI_Binary'))+
  scale_linetype_manual("",values=c('RLI or PDI' = 1,'randomized PDI' = 2))+
  labs(x = 'Year', y = 'Index value')+
  theme_bw()+
  theme(panel.grid =element_blank())



## combined PDI and randomised PDI of Mam
Plot_PDI_Mam.r <-ggplot()+
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_RLI.mean_noEX,colour = "RLI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_RLI.mean_noEX,colour = "RLI"),size=1)+
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_Issac.mean_noEX,colour = 'PDI_Issac',linetype="RLI or PDI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_Issac.mean_noEX,colour = 'PDI_Issac'),size=1)+
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_IUCN50.mean_noEX,colour = 'PDI_IUCN50',linetype="RLI or PDI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_IUCN50.mean_noEX,colour = 'PDI_IUCN50'),size=1)+
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_IUCN100.mean_noEX,colour = 'PDI_IUCN100',linetype="RLI or PDI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_IUCN100.mean_noEX,colour = 'PDI_IUCN100'),size=1)+
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_IUCN500.mean_noEX,colour = 'PDI_IUCN500',linetype="RLI or PDI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_IUCN500.mean_noEX,colour = 'PDI_IUCN500'),size=1)+
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_Pess.mean_noEX,colour = 'PDI_Pessimistic',linetype="RLI or PDI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_Pess.mean_noEX,colour = 'PDI_Pessimistic'),size=1)+
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_Binary.mean_noEX,colour = 'PDI_Binary',linetype="RLI or PDI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_Binary.mean_noEX,colour = 'PDI_Binary'),size=1)+
  
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_Issac.r.mean_noEX,colour = 'PDI_Issac.r',linetype="randomized PDI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_Issac.r.mean_noEX,colour = 'PDI_Issac.r'),size=1)+
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_IUCN50.r.mean_noEX,colour = 'PDI_IUCN50.r',linetype="randomized PDI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_IUCN50.r.mean_noEX,colour = 'PDI_IUCN50.r'),size=1)+
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_IUCN100.r.mean_noEX,colour = 'PDI_IUCN100.r',linetype="randomized PDI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_IUCN100.r.mean_noEX,colour = 'PDI_IUCN100.r'),size=1)+
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_IUCN500.r.mean_noEX,colour = 'PDI_IUCN500.r',linetype="randomized PDI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_IUCN500.r.mean_noEX,colour = 'PDI_IUCN500.r'),size=1)+
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_Pess.r.mean_noEX,colour = 'PDI_Pessimistic.r',linetype="randomized PDI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_Pess.r.mean_noEX,colour = 'PDI_Pessimistic.r'),size=1)+
  geom_line(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_Binary.r.mean_noEX,colour = 'PDI_Binary.r',linetype="randomized PDI"),size=0.5)+
  geom_point(data =  dd.m.c_noEX,aes(x = year.m,y = Mam_PDI_Binary.r.mean_noEX,colour = 'PDI_Binary.r'),size=1)+
  
  ylim(0.2,1.1)+scale_x_continuous(breaks=seq(1996, 2008, 4))+
  scale_colour_manual('', values = c('RLI' = "Red",'PDI_Issac' = 'pink','PDI_IUCN50' = 'steelblue1','PDI_IUCN100' = 'steelblue3',
                                     'PDI_IUCN500' = 'steelblue4','PDI_Pessimistic' = 'plum3','PDI_Binary' = 'plum4','PDI_Issac.r' = 'pink',
                                     'PDI_IUCN50.r' = 'steelblue1','PDI_IUCN100.r' = 'steelblue3',
                                     'PDI_IUCN500.r' = 'steelblue4','PDI_Pessimistic.r' = 'plum3','PDI_Binary.r' = 'plum4'),
                      breaks = c('RLI','PDI_Issac','PDI_IUCN50','PDI_IUCN100','PDI_IUCN500','PDI_Pessimistic','PDI_Binary'))+
  scale_linetype_manual("",values=c('RLI or PDI' = 1,'randomized PDI' = 2))+
  labs(x = 'Year', y = 'Index value')+
  theme_bw()+
  theme(panel.grid =element_blank())

pdf('E:/you know files/IC/PD project/Data2/fig_PDI_all.pdf', width = 8, height = 6)  
ggarrange(Plot_PDI_Coral, Plot_PDI_Mam,Plot_PDI_Coral.r, Plot_PDI_Mam.r, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)
dev.off()

pdf('E:/you know files/IC/PD project/Data2/fig_PDI_all2.pdf', width = 8, height = 6)  
ggarrange(Plot_PDI_Coral,Plot_PDI_Coral.r, Plot_PDI_Mam, Plot_PDI_Mam.r, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)
dev.off()

pdf('E:/you know files/IC/PD project/Data2/fig_PDI_all3.pdf', width = 8, height = 6)  
ggarrange(Plot_PDI_Coral.r,Plot_PDI_Mam.r, Plot_PDI_Coral,Plot_PDI_Mam, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)
dev.off()






