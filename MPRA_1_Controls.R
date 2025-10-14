## Code for Massively parallel functional screen identifies thousands of regulatory differences in human versus chimpanzee postcranial skeletal development

## Part 1: Library characteristics and controls


#load metadata
cartilage_metadata_experimental <- read.table(file = "~/Desktop/Capellini_Lab/MPRA/Cartilage_MPRA/cartilage_metadata_experimental_alignable.txt", header = T)

#get numbers and percentages of active controls

#create table for variables 
emvars_df <- data.frame(celltype = c("CHON002", "K562", "TC28"),
                        out_file = c("290k_CHON002_20250311.out", "290k_K562_20250311.out", "290k_TC28_20250311.out"),
                        glm_file = c("290k_CHON002_emVAR_glm_20250311.out", "290k_K562_emVAR_glm_20250311.out", "290k_TC28_emVAR_glm_20250311.out"),
                        active_controls = NA,
                        neg_controls = NA)

#iterate over cell type results data and extract control information
for(i in 1:3){
  temp_out_df <- read.table(emvars_df$out_file[i], header = T)
  #get number of positive control pairs and calculate percent success
  emvars_df$active_controls[i] <- temp_out_df %>% 
    dplyr::filter(grepl("emVarCtrl", project)) %>% 
    dplyr::filter(padj < 0.01) %>% 
    dplyr::filter(abs(log2FoldChange) > 1) %>% 
    dplyr::select(SNP) %>% unique() %>% nrow()/66*100
  #get number of negative control pairs and calculate percent success
  emvars_df$neg_controls[i] <- temp_out_df %>% 
    dplyr::filter(grepl("negCtrl", project)) %>% 
    dplyr::filter(padj < 0.01) %>% 
    dplyr::filter(abs(log2FoldChange) > 1) %>% 
    dplyr::select(SNP) %>% unique() %>% nrow()/71*100

}
#cleanup results
emvars_df$active_controls <- round(emvars_df$active_controls, digits = 1)

#create a figure showing activity/ 3 cell lines 
#based on Uebbing et al., PNAS, Fig.2 

attributesData <- read.table(file = "~/Desktop/290k_MPRA_Cluster_Results/AO_Cartilage_attributes.txt", header = T)
CHON002_normalized_counts <- read.table("~/Desktop/290k_MPRA_Cluster_Results/results/290k_20250311_CHON002_normalized_counts.out", header = T)
CHON002_normalized_counts_filtered <- CHON002_normalized_counts %>% 
  log2() %>% 
  rownames_to_column('ID') %>% 
  dplyr::filter(rownames(CHON002_normalized_counts) %in% attributesData$ID[which(attributesData$project == "SNP")]) %>% 
  dplyr::select(ID, c(CHON002_r1:CHON002_r5, plasmid_r1:plasmid_r4)) %>% 
  rowwise() %>%
  mutate(Min = min(c_across(CHON002_r1:plasmid_r4))) %>% 
  dplyr::filter(Min != -Inf) %>% 
  mutate(cDNA = median(CHON002_r1, CHON002_r2, CHON002_r3, CHON002_r4, CHON002_r5), 
         pDNA = median(plasmid_r1, plasmid_r2, plasmid_r3, plasmid_r4)) %>% 
  dplyr::select(ID, cDNA, pDNA)

K562_normalized_counts <- read.table("~/Desktop/290k_MPRA_Cluster_Results/results/290k_20250311_K562_normalized_counts.out", header = T)
K562_normalized_counts_filtered <- K562_normalized_counts %>% 
  log2() %>% 
  rownames_to_column('ID') %>% 
  dplyr::filter(rownames(K562_normalized_counts) %in% attributesData$ID[which(attributesData$project == "SNP")]) %>% 
  dplyr::select(ID, c(K562_r1:K562_r4, plasmid_r1:plasmid_r4)) %>% 
  rowwise() %>%
  mutate(Min = min(c_across(K562_r1:plasmid_r4))) %>% 
  dplyr::filter(Min != -Inf) %>% 
  mutate(cDNA = median(K562_r1, K562_r2, K562_r3, K562_r4), 
         pDNA = median(plasmid_r1, plasmid_r2, plasmid_r3, plasmid_r4)) %>% 
  dplyr::select(ID, cDNA, pDNA)

TC28_normalized_counts <- read.table("~/Desktop/290k_MPRA_Cluster_Results/results/290k_20250311_TC28_normalized_counts.out", header = T)
TC28_normalized_counts_filtered <- TC28_normalized_counts %>% 
  log2() %>% 
  rownames_to_column('ID') %>% 
  dplyr::filter(rownames(TC28_normalized_counts) %in% attributesData$ID[which(attributesData$project == "SNP")]) %>% 
  dplyr::select(ID, c(TC28_r1:TC28_r5, plasmid_r1:plasmid_r4)) %>% 
  rowwise() %>%
  mutate(Min = min(c_across(TC28_r1:plasmid_r4))) %>% 
  dplyr::filter(Min != -Inf) %>% 
  mutate(cDNA = median(TC28_r1, TC28_r2, TC28_r3, TC28_r4, TC28_r5), 
         pDNA = median(plasmid_r1, plasmid_r2, plasmid_r3, plasmid_r4)) %>% 
  dplyr::select(ID, cDNA, pDNA)

#add metadata
add_norm_count_meta_columns <- function(filtered_counts, glm_file, celltype){
  filtered_counts$Species <- "Human"
  filtered_counts$Species[grep(pattern = "_chimp", x = filtered_counts$ID)] <- "Chimp"
  filtered_counts$Celltype <- celltype
  filtered_counts$Status <- "Inactive"
  #load glm results to identify emVars
  temp_glm_df <- read.table(glm_file, header = T)
  
  #get counts of paired oligos
  human_active_pairs <- temp_glm_df %>% 
    dplyr::filter((A_logPadj_BH > 2 & abs(A_log2FC) > 1 )) %>% pull(SNP)
  human_active_pairs2 <- paste(human_active_pairs, "_hs", sep = "")
  chimp_active_pairs <- temp_glm_df %>% 
    dplyr::filter((B_logPadj_BH > 2 & abs(B_log2FC) > 1 )) %>% pull(SNP)
  chimp_active_pairs2 <- paste(chimp_active_pairs, "_chimp", sep = "")
  #get counts of emVars
  human_emVars <- temp_glm_df %>% 
    dplyr::filter((A_logPadj_BH > 2 & abs(A_log2FC) > 1 )) %>% 
    #the first filter is unnecessary but makes me feel better
    dplyr::filter(Skew_logFDR_act > 1 & Log2Skew < 0) %>% 
    pull(SNP)
  human_emVars2 <- paste(human_emVars, "_hs", sep = "")
  chimp_emVars <- temp_glm_df %>% 
    dplyr::filter((B_logPadj_BH > 2 & abs(B_log2FC) > 1 )) %>% 
    #the first filter is unnecessary but makes me feel better
    dplyr::filter(Skew_logFDR_act > 1 & Log2Skew > 0) %>% 
    pull(SNP)
  chimp_emVars2 <- paste(chimp_emVars, "_chimp", sep = "")
  
  filtered_counts$Status[which(filtered_counts$ID %in% human_active_pairs2)] <- "Human-active"
  filtered_counts$Status[which(filtered_counts$ID %in% chimp_active_pairs2)] <- "Chimp-active"
  filtered_counts$Status[which(filtered_counts$ID %in% human_emVars2)] <- "Human-biased"
  filtered_counts$Status[which(filtered_counts$ID %in% chimp_emVars2)] <- "Chimp-biased"
  return(filtered_counts)
}

#add metadata to each normalized count file and combine
CHON002_normalized_counts_filtered2 <- add_norm_count_meta_columns(filtered_counts = CHON002_normalized_counts_filtered, 
                                                                   celltype = "CHON002", 
                                                                   glm_file = "290k_CHON002_emVAR_glm_20250311.out")

K562_normalized_counts_filtered2 <- add_norm_count_meta_columns(filtered_counts = K562_normalized_counts_filtered, 
                                                                celltype = "K562", 
                                                                glm_file = "290k_K562_emVAR_glm_20250311.out")

TC28_normalized_counts_filtered2 <- add_norm_count_meta_columns(filtered_counts = TC28_normalized_counts_filtered, 
                                                                celltype = "TC28", 
                                                                glm_file = "290k_TC28_emVAR_glm_20250311.out")

#bind by rows
normalized_counts_filtered2 <- rbind(CHON002_normalized_counts_filtered2, K562_normalized_counts_filtered2, TC28_normalized_counts_filtered2)

normalized_counts_filtered2$Species <- factor(normalized_counts_filtered2$Species, levels = c("Human", "Chimp"))
counts_plot <- ggplot(data = normalized_counts_filtered2, aes(x = pDNA, y = cDNA, color = Status)) + 
  geom_point(size = 0.25) + 
  facet_grid(rows =vars(Species), cols = vars(Celltype)) + 
  scale_color_manual(, values = c("#44AA99", "#117733","#88CCEE","#0079B2", "#888888"))+
  labs(x = "Median log2 pDNA counts", y = "Median log2 cDNA counts") + 
  guides(colour = guide_legend(override.aes = list(size=10))) +
  theme_bw() + 
  theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#save results
ggsave(plot = counts_plot, 
       filename = "~/Dropbox/Cartilage MPRA Paper/Code/S2_counts_active_emvar_plot.png", 
       device = "png", dpi = 300, height = 10, width = 10, units = "in")