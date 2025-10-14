#build attributes file for 290k Cartilage MPRA library

library(stringr)
countsData <- read.table("~/Desktop/290k_MPRA_Cluster_Results/AO_Cartilage.count", header = T)

#create attributes dataframe
attributes <- data.frame(id = unique(countsData$Oligo), 
                         snp = unique(countsData$Oligo), 
                         chromosome = NA, 
                         snp_pos = NA, 
                         ref_allele = NA, 
                         alt_allele = NA, 
                         allele = NA, 
                         window = "center",  
                         strand ="fwd",  
                         project = "SNP")

#load controls for identification
#read in controls data
MPRA_270bp_pos_controls <- read.table("~/Desktop/Capellini_Lab/MPRA/Cartilage_MPRA/MPRA_controls_active.id.270.seq")
names(MPRA_270bp_pos_controls) <- c("ID","sequence")

MPRA_270bp_neg_controls <- read.table("~/Desktop/Capellini_Lab/MPRA/Cartilage_MPRA/K562_AllTREs_110819_Yu_ORF-NegCtrl_270bp.seq")
names(MPRA_270bp_neg_controls) <- c("ID","sequence")

#add control data
attributes$project[which(attributes$id %in% MPRA_270bp_pos_controls$ID)] <- "emVarCtrl"
attributes$project[which(attributes$id %in% MPRA_270bp_neg_controls$ID)] <- "negCtrl"

#add positional information
split1 <- as.data.frame(unlist(str_split_fixed(attributes$id[which(attributes$project == "SNP")], pattern = ":", n = 2)))
attributes$chromosome[which(attributes$project == "SNP")] <- split1$V1
split2 <- as.data.frame(unlist(str_split_fixed(split1$V2, pattern = "_", n = 3)))

#add SNP_ID
attributes$snp[which(attributes$project == "SNP")] <- paste(split1$V1, ":", split2$V1, "_", split2$V2, sep = "")
#correct start position of second tiles
split2$V1 <- as.numeric(split2$V1)
split2$V1[which(split2$V2 == "2")] <- split2$V1[which(split2$V2 == "2")] + 269
attributes$snp_pos[which(attributes$project == "SNP")] <- split2$V1
#add window information
attributes$window[which(attributes$project == "SNP")][which(split2$V2 == "1")] <- "1"
attributes$window[which(attributes$project == "SNP")][which(split2$V2 == "2")] <- "2"

#add allele information

attributes$allele[grep(x = attributes$id, pattern = "_hs")] <- "ref"
attributes$allele[grep(x = attributes$id, pattern = "_chimp")] <- "alt"

#I think there is a capitalization typo...
attributes$ID <- attributes$id
attributes$SNP <- attributes$snp
attributes$chr <- attributes$chromosome

write.table(x = attributes, file = "~/Desktop/290k_MPRA_Cluster_Results/AO_Cartilage_attributes.txt", col.names = T)

