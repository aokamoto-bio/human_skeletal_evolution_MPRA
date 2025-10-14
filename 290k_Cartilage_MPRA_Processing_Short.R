setwd("~/Desktop/290k_MPRA_Cluster_Results/")
#load packages
library("MPRAmodel")
library("tidyverse")

countsData <- read.table("~/Desktop/290k_MPRA_Cluster_Results/AO_Cartilage.count", header = T)
attributesData <- read.table(file = "~/Desktop/290k_MPRA_Cluster_Results/AO_Cartilage_attributes.txt", header = T)

#remove problematic replicate
countsData$K562_r5 <- NULL

conditionData <- read.table("~/Desktop/290k_MPRA_Cluster_Results/AO_Cartilage_condition.txt", header = F, row.names = 1)
conditionData <- conditionData %>% slice(-10)

filePrefix <- "290k"
projectName <- "290k_MPRA"
exclList <- c()
plotSave <- T
altRef <- T 
method <- 'ss'
tTest <- F
DEase <- T
cSkew <- T 
correction <- "BH"
cutoff <- 0.01
emVAR_cutoff <- 0.1
upDisp <- T
prior <- F
raw <- T
paired <- F
negCtrlName <- "negCtrl"

# Make sure that the plots and results directories are present in the current directory
mainDir <- getwd()
#dir.create(file.path(mainDir, "plots"), showWarnings = FALSE)
#dir.create(file.path(mainDir, "results"), showWarnings = FALSE)
# Resolve any multi-project conflicts, run normalization, and write celltype specific results files
attributesData <- addHaplo(attributesData, negCtrlName, posCtrlName, projectName)
message("running DESeq")
analysis_out <- dataOut(countsData, 
                        attributesData, 
                        conditionData, 
                        altRef = altRef, 
                        exclList = exclList, 
                        file_prefix = filePrefix, 
                        method = method, 
                        negCtrlName = "negCtrl", 
                        tTest = tTest, 
                        DEase = T, 
                        cSkew = T, 
                        correction = correction, 
                        cutoff = 0.01, 
                        upDisp = T, 
                        prior = F, 
                        paired = F)
