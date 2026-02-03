## Code for Massively parallel functional screen identifies thousands of regulatory differences in human versus chimpanzee postcranial skeletal development

## Part 2: Characteristics of active sequences

#load data
#modify to specify file path for your computer
setwd("~/Desktop/290k_MPRA_Cluster_Results/results/")
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggVennDiagram)
library(cowplot)
library(magick) # for adding images to ggplots
#load custom functions for analysis
#this code is available from https://github.com/aokamoto-bio/Human_Autopod_Evolution
source("~/Desktop/Capellini_Lab/Weekly_Coding/RNA_Analysis_Functions.R") 

#load metadata
cartilage_metadata_experimental <- read.table(file = "~/Desktop/Capellini_Lab/MPRA/Cartilage_MPRA/cartilage_metadata_experimental_alignable.txt", header = T)

#load active pairs for each cell type
CHON002_active <- read.table("290k_CHON002_emVAR_glm_20250311.out", header = T) %>% 
  dplyr::filter((A_logPadj_BH > 2 & abs(A_log2FC) > 1 ) | 
                  (B_logPadj_BH > 2 & abs(B_log2FC) > 1 )) %>% 
  merge(y = cartilage_metadata_experimental, by = "SNP", all.x = F)

TC28_active <- read.table("290k_TC28_emVAR_glm_20250311.out", header = T) %>% 
  dplyr::filter((A_logPadj_BH > 2 & abs(A_log2FC) > 1 ) | 
                  (B_logPadj_BH > 2 & abs(B_log2FC) > 1 )) %>% 
  merge(y = cartilage_metadata_experimental, by = "SNP", all.x = F, all.y = F)

K562_active <- read.table("290k_K562_emVAR_glm_20250311.out", header = T) %>% 
  dplyr::filter((A_logPadj_BH > 2 & abs(A_log2FC) > 1 ) | 
                  (B_logPadj_BH > 2 & abs(B_log2FC) > 1 )) %>% 
  merge(y = cartilage_metadata_experimental, by = "SNP", all.x = F)

#load emVar pairs for each cell type
CHON002_emVars <- read.table("290k_CHON002_emVAR_glm_20250311.out", header = T) %>% 
  dplyr::filter((A_logPadj_BH > 2 & abs(A_log2FC) > 1 ) | 
                  (B_logPadj_BH > 2 & abs(B_log2FC) > 1 )) %>% 
  dplyr::filter(Skew_logFDR_act > 1) %>% 
  merge(y = cartilage_metadata_experimental, by = "SNP", all.x = F)

TC28_emVars <- read.table("290k_TC28_emVAR_glm_20250311.out", header = T) %>% 
  dplyr::filter((A_logPadj_BH > 2 & abs(A_log2FC) > 1 ) | 
                  (B_logPadj_BH > 2 & abs(B_log2FC) > 1 )) %>% 
  dplyr::filter(Skew_logFDR_act > 1) %>% 
  merge(y = cartilage_metadata_experimental, by = "SNP", all.x = F)

K562_emVars <- read.table("290k_K562_emVAR_glm_20250311.out", header = T) %>% 
  dplyr::filter((A_logPadj_BH > 2 & abs(A_log2FC) > 1 ) | 
                  (B_logPadj_BH > 2 & abs(B_log2FC) > 1 )) %>% 
  dplyr::filter(Skew_logFDR_act > 1) %>% 
  merge(y = cartilage_metadata_experimental, by = "SNP", all.x = F)

#count total sequence pairs analyzed across dataset
(total_counts <- read.table("~/Desktop/290k_MPRA_Cluster_Results/results/290k_CHON002_emVAR_glm_20250311.out", header = T) %>% 
  dplyr::filter(SNP %in% cartilage_metadata_experimental$SNP) %>% 
  nrow())

#count total regions analyzed across dataset
(total_regions_counts <- read.table("~/Desktop/290k_MPRA_Cluster_Results/results/290k_CHON002_emVAR_glm_20250311.out", header = T) %>% 
  dplyr::filter(SNP %in% cartilage_metadata_experimental$SNP) %>% 
  merge(y = cartilage_metadata_experimental, by = "SNP", all.x = F) %>% 
  pull(region) %>% unique() %>% length())

#total regions with activity in at least one cell type
(active_regions_count <- length(unique(c(CHON002_active$region, K562_active$region, TC28_active$region))))
#percent
active_regions_count/total_regions_counts*100

#total active regions per cell type
(CHON_active_regions_count <- length(unique(CHON002_active$region)))
(TC28_active_regions_count <- length(unique(TC28_active$region)))
(K562_active_regions_count <- length(unique(K562_active$region)))

#plot venn diagram to visualize active region sharing between cell types
activity_venn <- ggVennDiagram(
  x = list(unique(CHON002_active$region), unique(K562_active$region), unique(TC28_active$region)),
  category.names = c("        CHON002" , "K562" , "TC28"), 
  label_alpha = 0, label_size = 3
) + 
  labs(fill = "Count") + 
  theme(plot.title=element_text(hjust=0.5)) + 
  scale_fill_gradient(low="grey90",high = "red")

### activity sharing between tiles section

#count total number of regions with both tiles active
(all_2tile_count <- rbind(CHON002_active, K562_active, TC28_active) %>% 
  dplyr::select(region, tile) %>% 
  unique() %>% 
  group_by(region) %>% 
  dplyr::summarise(n =n()) %>% 
  dplyr::filter(n >1) %>% nrow())

#count total number of regions with only one tile active
(all_1tile_count <- rbind(CHON002_active, K562_active,TC28_active) %>% 
  dplyr::select(region, tile) %>% 
  unique() %>% 
  group_by(region) %>% 
  summarise(n =n()) %>% 
  dplyr::filter(n == 1) %>% nrow())

#dive deeper into two tile regions
all_2tile <- rbind(CHON002_active, K562_active, TC28_active) %>% 
  dplyr::select(region, tile) %>% 
  unique() %>% 
  group_by(region) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n > 1) %>% pull(region)

chon_2tile <- rbind(CHON002_active) %>% 
  dplyr::select(region, tile) %>% 
  unique() %>% 
  group_by(region) %>% 
  dplyr::summarise(n =n()) %>% 
  dplyr::filter(n > 1) %>% pull(region)

k562_2tile <- rbind(K562_active) %>% 
  dplyr::select(region, tile) %>% 
  unique() %>% 
  group_by(region) %>% 
  dplyr::summarise(n =n()) %>% 
  dplyr::filter(n > 1) %>% pull(region)

tc28_2tile <- rbind(TC28_active) %>% 
  dplyr::select(region, tile) %>% 
  unique() %>% 
  group_by(region) %>% 
  dplyr::summarise(n =n()) %>% 
  dplyr::filter(n > 1) %>% pull(region)

#count regions with 2 regions across two cell types
length(all_2tile) - length(unique(c(tc28_2tile, chon_2tile, k562_2tile)))
tiles_split <- all_2tile[!all_2tile %in% unique(c(tc28_2tile, chon_2tile, k562_2tile))]
length(intersect(tiles_split[which(tiles_split %in% TC28_active$region)], tiles_split[which(tiles_split %in% CHON002_active$region)]))
length(intersect(tiles_split[which(tiles_split %in% TC28_active$region)], tiles_split[which(tiles_split %in% K562_active$region)]))
length(intersect(tiles_split[which(tiles_split %in% CHON002_active$region)], tiles_split[which(tiles_split %in% K562_active$region)]))

TC28_active_tiles <- TC28_active %>% 
  dplyr::filter(region %in% tc28_2tile) %>% 
  dplyr::select(region, tile, A_log2FC, B_log2FC) %>% distinct() %>% 
  group_by(region) %>% 
  pivot_wider(values_from = c(A_log2FC, B_log2FC), names_from = tile)

CHON002_active_tiles <- CHON002_active %>% 
  dplyr::filter(region %in% chon_2tile) %>% 
  dplyr::select(region, tile, A_log2FC, B_log2FC) %>% distinct() %>% 
  group_by(region) %>% 
  pivot_wider(values_from = c(A_log2FC, B_log2FC), names_from = tile)

K562_active_tiles <- K562_active %>% 
  dplyr::filter(region %in% k562_2tile) %>% 
  dplyr::select(region, tile, A_log2FC, B_log2FC) %>% distinct() %>% 
  group_by(region) %>% 
  pivot_wider(values_from = c(A_log2FC, B_log2FC), names_from = tile)

#build scatterplots to show results
two_tile_simple_corr_hum <- data.frame(tile1_skew = c(TC28_active_tiles$A_log2FC_1, CHON002_active_tiles$A_log2FC_1, K562_active_tiles$A_log2FC_1), 
                                   tile2_skew = c(TC28_active_tiles$A_log2FC_2, CHON002_active_tiles$A_log2FC_2, K562_active_tiles$A_log2FC_2), 
                                   species = "Human", 
                                   cellline = c(rep("TC28", times = nrow(TC28_active_tiles)), 
                                                rep("CHON002", times = nrow(CHON002_active_tiles)), 
                                                rep("K562", times = nrow(K562_active_tiles)))
                                   )

two_tile_simple_corr_chimp <- data.frame(tile1_skew = c(TC28_active_tiles$B_log2FC_1, CHON002_active_tiles$B_log2FC_1, K562_active_tiles$B_log2FC_1), 
                                       tile2_skew = c(TC28_active_tiles$B_log2FC_2, CHON002_active_tiles$B_log2FC_2, K562_active_tiles$B_log2FC_2), 
                                       species = "Chimp", 
                                       cellline = c(rep("TC28", times = nrow(TC28_active_tiles)), 
                                                    rep("CHON002", times = nrow(CHON002_active_tiles)), 
                                                    rep("K562", times = nrow(K562_active_tiles)))
                                    )
#combine species 
two_tile_simple_corr <- rbind(two_tile_simple_corr_hum, two_tile_simple_corr_chimp)
#factor to list human first
two_tile_simple_corr$species <- factor(two_tile_simple_corr$species, levels = c("Human", "Chimp"))

#create the plot
two_tile_simple_corr_plot <- ggplot(data = two_tile_simple_corr, 
       mapping = aes(x = tile1_skew, y = tile2_skew)) + 
  geom_point(size = 0.2) + 
  facet_grid(species~cellline) + 
  stat_cor(method="pearson", label.y = 12.5) + 
  labs(x = "Tile 1 log2 fold change activity", y = "Tile 2 log2 fold change activity") +
  ylim(c(-20, 15))

#now for the fancier case
#compare pairwise sets 

split_tiles_cell_type_active <- function(tiles_split_list, cell_df1, cell_df2){
  
  #get list of split tiles in these two dfs
  split_tiles <- intersect(tiles_split_list[which(tiles_split_list %in% cell_df1$region)], tiles_split_list[which(tiles_split_list %in% cell_df2$region)])
  
  split_tiles_df <- data.frame(SNP = cartilage_metadata_experimental$SNP[which(cartilage_metadata_experimental$region %in% split_tiles)],
                               region = cartilage_metadata_experimental$region[which(cartilage_metadata_experimental$region %in% split_tiles)], 
                               tile = cartilage_metadata_experimental$tile[which(cartilage_metadata_experimental$region %in% split_tiles)], 
                               A_log2FC = NA, 
                               B_log2FC = NA)
  
  #add emvar tile info from relevant cell type
  split_tiles_df$A_log2FC[which(split_tiles_df$SNP %in% cell_df1$SNP)] <- cell_df1$A_log2FC[which(cell_df1$SNP %in% split_tiles_df$SNP)]
  split_tiles_df$B_log2FC[which(split_tiles_df$SNP %in% cell_df1$SNP)] <- cell_df1$B_log2FC[which(cell_df1$SNP %in% split_tiles_df$SNP)]
  split_tiles_df$A_log2FC[which(split_tiles_df$SNP %in% cell_df2$SNP)] <- cell_df2$A_log2FC[which(cell_df2$SNP %in% split_tiles_df$SNP)]
  split_tiles_df$B_log2FC[which(split_tiles_df$SNP %in% cell_df2$SNP)] <- cell_df2$B_log2FC[which(cell_df2$SNP %in% split_tiles_df$SNP)]
  
  split_tiles_df2 <- split_tiles_df %>% 
    dplyr::select(-SNP) %>% 
    dplyr::filter(region %in% split_tiles) %>% 
    dplyr::select(region, tile, A_log2FC, B_log2FC) %>% 
    distinct() %>% 
    group_by(region) %>% 
    pivot_wider(values_from = c(A_log2FC, B_log2FC), names_from = tile) %>% 
    drop_na()
  
  #return filled in dataframe
  return(split_tiles_df2)
}

#generate split tiles dfs

#chon versus tc28
split_tiles_df_chon_tc28 <- split_tiles_cell_type_active(tiles_split_list = tiles_split, cell_df1 = TC28_active, cell_df2 = CHON002_active)
#chon versus k562
split_tiles_df_chon_k562 <- split_tiles_cell_type_active(tiles_split_list = tiles_split, cell_df1 = K562_active, cell_df2 = CHON002_active)
#tc28 versus k562
split_tiles_df_tc28_k562 <- split_tiles_cell_type_active(tiles_split_list = tiles_split, cell_df1 = K562_active, cell_df2 = TC28_active)

#build scatterplots to show results
two_tile_complex_corr_hum <- data.frame(tile1_skew = c(split_tiles_df_chon_tc28$A_log2FC_1, split_tiles_df_chon_k562$A_log2FC_1, split_tiles_df_tc28_k562$A_log2FC_1), 
                                       tile2_skew = c(split_tiles_df_chon_tc28$A_log2FC_2, split_tiles_df_chon_k562$A_log2FC_2, split_tiles_df_tc28_k562$A_log2FC_2), 
                                       species = "Human", 
                                       cellline = c(rep("CHON002 vs TC28", times = nrow(split_tiles_df_chon_tc28)), 
                                                    rep("CHON002 vs K562", times = nrow(split_tiles_df_chon_k562)), 
                                                    rep("TC28 vs K562", times = nrow(split_tiles_df_tc28_k562)))
)

two_tile_complex_corr_chimp <- data.frame(tile1_skew = c(split_tiles_df_chon_tc28$B_log2FC_1, split_tiles_df_chon_k562$B_log2FC_1, split_tiles_df_tc28_k562$B_log2FC_1), 
                                         tile2_skew = c(split_tiles_df_chon_tc28$B_log2FC_2, split_tiles_df_chon_k562$B_log2FC_2, split_tiles_df_tc28_k562$B_log2FC_2), 
                                         species = "Chimp", 
                                         cellline = c(rep("CHON002 vs TC28", times = nrow(split_tiles_df_chon_tc28)), 
                                                      rep("CHON002 vs K562", times = nrow(split_tiles_df_chon_k562)), 
                                                      rep("TC28 vs K562", times = nrow(split_tiles_df_tc28_k562)))
)

#combine species 
two_tile_complex_corr <- rbind(two_tile_complex_corr_hum, two_tile_complex_corr_chimp)
#factor to list human first
two_tile_complex_corr$species <- factor(two_tile_complex_corr$species, levels = c("Human", "Chimp"))

#create the plot
two_tile_complex_corr_plot <- ggplot(data = two_tile_complex_corr, 
                                    mapping = aes(x = tile1_skew, y = tile2_skew)) + 
  geom_point(size = 0.2) + 
  facet_grid(species~cellline) + 
  stat_cor(method="pearson", label.y = 12) + 
  labs(x = "Tile 1 log2 fold change activity", y = "Tile 2 log2 fold change activity") +
  ylim(c(-20, 15))

#create the supplementary figure

two_tile_figure <- ggdraw() +
  draw_plot(two_tile_simple_corr_plot, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(two_tile_complex_corr_plot, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B"), 
                  size = 12,
                  x = c(0, 0), 
                  y = c(1, 0.5)) + 
  bgcolor("white")

ggsave(plot = two_tile_figure , 
       filename = "~/Dropbox/Cartilage MPRA Paper/Code/S3_two_tile_scatterplots.png", 
       device = "png", dpi = 300, height = 8, width = 6, units = "in")

#get split tiles 
#may not be used? 
TC28_split_tiles <- TC28_active %>% 
  dplyr::filter(region %in% tiles_split) %>% 
  dplyr::select(region, tile, A_log2FC, B_log2FC) %>% distinct() %>% 
  group_by(region) %>% 
  pivot_wider(values_from = c(A_log2FC, B_log2FC), names_from = tile)

CHON002_split_tiles <- CHON002_active %>% 
  dplyr::filter(region %in% tiles_split) %>% 
  dplyr::select(region, tile, A_log2FC, B_log2FC) %>% distinct() %>% 
  group_by(region) %>% 
  pivot_wider(values_from = c(A_log2FC, B_log2FC), names_from = tile)

K562_split_tiles <- K562_active %>% 
  dplyr::filter(region %in% tiles_split) %>% 
  dplyr::select(region, tile, A_log2FC, B_log2FC) %>% distinct() %>% 
  group_by(region) %>% 
  pivot_wider(values_from = c(A_log2FC, B_log2FC), names_from = tile)

#plot pie chart of active sequences

#format data
activity_pie_data <- data.frame(
  class = c("Inactive", "Active (2 tile)", "Active (1 tile)"), 
  count = c(total_regions_counts - active_regions_count, all_2tile_count, all_1tile_count)
)

#format for plotting
activity_pie_data2 <- activity_pie_data %>% 
  mutate(pct = count/sum(count)*100)

activity_pie_data2$label = c("37,250\n(54.8%)","4,529\n(6.6%)",  "26,207\n(38.5%)")

#generate plot
activity_pie <- ggplot(data = activity_pie_data2, aes(x="", y = pct, fill = class)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  #geom_label_repel(aes(label = count, fill = NULL),
  #                 size = 1.5, nudge_x = 1, show.legend = FALSE)+
  geom_text(aes(label = label, x = 1.1), position = position_stack(vjust=0.5), size = 3) +
  theme_void() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#88CCEE", "lightcyan1","grey"), name = NULL) + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

#test for bias towards human or chimp activity
binom.test(length(which(CHON002_active$A_logPadj_BH > 2 )), n = (length(which(CHON002_active$B_logPadj_BH > 2)) + length(which(CHON002_active$A_logPadj_BH > 2 ))), p = 0.5)
binom.test(length(which(TC28_active$A_logPadj_BH > 2 )), n = (length(which(TC28_active$B_logPadj_BH > 2)) + length(which(TC28_active$A_logPadj_BH > 2 ))), p = 0.5)
binom.test(length(which(K562_active$A_logPadj_BH > 2 )), n = (length(which(K562_active$B_logPadj_BH > 2)) + length(which(K562_active$A_logPadj_BH > 2 ))), p = 0.5)

get_MPRA_bed <- function(region_IDs, filename, type = "region"){
  filename_full <- paste("~/Desktop/290k_MPRA_Cluster_Results/beds/", filename, "_hg38.bed", sep ="")
  if(toupper(type) == "REGION"){
    cartilage_metadata_experimental %>% 
      dplyr::filter(region %in% region_IDs) %>% 
      mutate(region_stop = region_start + 540) %>% 
      dplyr::select(chr, region_start, region_stop) %>% 
      unique() %>% 
      write_bed(filename = filename_full, ncol = 3)
  }
  if(toupper(type) == "TILE"){
    cartilage_metadata_experimental %>% 
      dplyr::filter(SNP %in% region_IDs) %>% 
      dplyr::select(chr, start, end, SNP) %>% 
      unique() %>% 
      write_bed(filename = filename_full, ncol = 4)
  }
  return(filename_full)
}

# create a function to visualize a bed file on the human skeleton
# Alexander Okamoto
# March 14, 2025

#define function
#function inputs
#bed_file: file for overlap
#ATAC_path: file path to (and including) the "Capellini Beds" directory containing the human ATAC files, results table, heatmap function, 
#default ATAC_path is for Alexander's computer
#timepoint: default is both (E54 & E67), accepts "E54"/"early" or "E67"/"late
#also accept "unique" which looks at unique regions identified by merging timepoints
#genome: default hg38, also accepts hg19
#label: label for resulting plot (and partial plot file name)

bed_capellini_skeleton_overlaps <- function(
    bed_file, 
    ATAC_path = "~/Desktop/Capellini_Lab/Weekly_Coding/Capellini_Beds/", 
    timepoint = "both", 
    genome = "hg38", 
    label) {
  
  #load packages
  require(tidyverse)
  
  #load dataframe to store overlaps and with bed files path
  results_file <- paste(ATAC_path, "capellini_skeleton_heatmap_beds.txt", sep ="")
  
  #check that file paths work
  if(!file.exists(results_file)){
    stop("SKELETAL BEDS NOT FOUND. PLEASE CHECK FILEPATH TO 'Capellini_Beds'!!!")
    
  }
  results_df <- read.table(results_file, header = T)
  
  #load skeletal heatmap function
  source(paste(ATAC_path, "Skeleton_Heatmap_V1.R", sep =""))
  
  #first, check desired genome and assign appropriate files
  if(toupper(genome) == "HG38"){
    print("Using hg38 genome coordinates")
    #add appropriate file path
    results_df$early_file <- paste(ATAC_path, results_df$early_file_hg38, sep = "")
    results_df$late_file <- paste(ATAC_path, results_df$late_file_hg38, sep = "")
    results_df$unique_file <- paste(ATAC_path, results_df$unique_hg38, sep = "")
    
  }
  
  if(toupper(genome) == "HG19"){
    print("Using hg19 genome coordinates")
    #add appropriate file path
    results_df$early_file <- paste(ATAC_path, results_df$early_file_hg19, sep = "")
    results_df$late_file <- paste(ATAC_path, results_df$late_file_hg19, sep = "")
    results_df$unique_file <- paste(ATAC_path, results_df$unique_hg19, sep = "")
  }
  
  #next, check if only one or both timepoints desired
  if(toupper(timepoint) == "BOTH"){
    print("Performing overlaps at both timepoints")
    #iterate over tissues
    for(x in 1:nrow(results_df)){
      if(!x %in% 29:30){
        results_df$value[x] <- system(paste("cat ", results_df$early_file[x], " ", results_df$late_file[x], " | cut -f'1-3' | sort -k 1,1 -k2,2n | bedtools merge -i stdin | bedtools intersect -a stdin -b ", bed_file, " | wc -l", sep =""), intern =T) %>% 
          readr::parse_number()
      }
      if(x %in% 29:30){
        results_df$value[x] <- system(paste("cut -f'1-3' ", results_df$late_file[x], " | bedtools intersect -a stdin -b ", bed_file, " | wc -l", sep =""), intern =T) %>% 
          readr::parse_number()
      }
    } #end iterate over tissues
  } #end both timepoint option
  
  if(toupper(timepoint) == "EARLY" | toupper(timepoint) == "E54"){
    print("Performing overlaps at both timepoints")
    #iterate over tissues
    for(x in 1:nrow(results_df)){
      if(!x %in% 29:30){
        results_df$value[x] <- system(paste("cut -f'1-3' ", results_df$early_file[x], " | bedtools intersect -a stdin -b ", bed_file, " | wc -l", sep =""), intern =T) %>% 
          readr::parse_number()
      }
    } #end iterate over tissues
  } #end early option
  
  if(toupper(timepoint) == "LATE" | toupper(timepoint) == "E67"){
    print("Performing overlaps at both timepoints")
    #iterate over tissues
    for(x in 1:nrow(results_df)){
      results_df$value[x] <- system(paste("cut -f'1-3' ", results_df$late_file[x], " | bedtools intersect -a stdin -b ", bed_file, " | wc -l", sep =""), intern =T) %>% 
        readr::parse_number()
    } #end iterate over tissues
  } #end early option
  
  if(toupper(timepoint) == "UNIQUE"){
    print("Performing overlaps for unique regions")
    #iterate over tissues
    for(x in 1:nrow(results_df)){
      results_df$value[x] <- system(paste("cut -f'1-3' ", results_df$unique_file[x], " | bedtools intersect -a stdin -b ", bed_file, " | wc -l", sep =""), intern =T) %>% 
        readr::parse_number()
    } #end iterate over tissues
  } #end early option
  
  try(system("mkdir ~/Desktop/Autopod_Skeleton_Figures"))
  print("results in folder: ~/Desktop/Autopod_Skeleton_Figures")
  plot_skeleton_heatmap(expr_df = results_df %>% dplyr::select(tissue, value), 
                        inc.autopods = T, 
                        file_name = paste("~/Desktop/Autopod_Skeleton_Figures/", label, "_", timepoint, "_skeleton_heatmap", sep = ""), 
                        key_label = label, 
                        binary = F)
  return(results_df %>% dplyr::select(tissue, value))
}



#create wrapper function to run MPRA bed function efficiently
MPRA_multi_bed_overlaps <- function(emvar_regions, active_regions, filename, label){
  #first get the total overlaps for plotting
  all_bed <- get_MPRA_bed(cartilage_metadata_experimental$region, filename = "290k_all_bed")
  all_tested_oligo_overlaps <- bed_capellini_skeleton_overlaps(bed_file = all_bed, label = "all")
  all_tested_oligo_unique_overlaps <- bed_capellini_skeleton_overlaps(bed_file = all_bed, label = "unique", timepoint = "unique")
  #modify column names for clarity and to allow merging
  colnames(all_tested_oligo_overlaps)[2] <- "total_overlaps"
  colnames(all_tested_oligo_unique_overlaps)[2] <- "unique_overlaps"
  
  #first for active regions
  active_bed <- get_MPRA_bed(active_regions, 
                             filename = paste(filename, "_active", sep =""))
  active_overlaps <- bed_capellini_skeleton_overlaps(bed_file = active_bed, 
                                                     label = paste(label, "_active", sep =""), 
                                                     timepoint = "both")
  active_overlaps_unique <- bed_capellini_skeleton_overlaps(bed_file = active_bed, 
                                                            label = paste(label, "_active_uniq", sep =""), 
                                                            timepoint = "unique")
  #second for emvars regions
  emvar_bed <- get_MPRA_bed(active_regions, 
                            filename = paste(filename, "_emVars", sep =""))
  emvar_overlaps <- bed_capellini_skeleton_overlaps(bed_file = emvar_bed, 
                                                    label = paste(label, "_emVars", sep =""), 
                                                    timepoint = "both")
  emvar_overlaps_unique <- bed_capellini_skeleton_overlaps(bed_file = emvar_bed, 
                                                           label = paste(label, "_emVars_uniq", sep =""), 
                                                           timepoint = "unique")
  #modify column names for clarity and to allow merging
  colnames(active_overlaps)[2] <- "active_overlaps"
  colnames(active_overlaps_unique)[2] <- "active_unique_overlaps"
  colnames(emvar_overlaps)[2] <- "emvar_overlaps"
  colnames(emvar_overlaps_unique)[2] <- "emvar_unique_overlaps"
  
  overlap_sharing <- merge(all_tested_oligo_overlaps, all_tested_oligo_unique_overlaps) %>% 
    merge(active_overlaps) %>% 
    merge(active_overlaps_unique)%>% 
    merge(emvar_overlaps) %>% 
    merge(emvar_overlaps_unique) 
  overlap_sharing$celltype  <- label
  return(overlap_sharing)
} #end function

all_overlap_sharing <- MPRA_multi_bed_overlaps(emvar_regions = unique(c(CHON002_emVars$region, TC28_emVars$region, K562_emVars$region)), 
                                               active_regions = unique(c(CHON002_active$region, TC28_active$region, K562_active$region)), 
                                               filename = "290k_ALL", 
                                               label = "ALL")

CHON_overlap_sharing <- MPRA_multi_bed_overlaps(emvar_regions = CHON002_emVars$region, 
                                                active_regions = CHON002_active$region, 
                                                filename = "290k_CHON", 
                                                label = "CHON")
TC28_overlap_sharing <- MPRA_multi_bed_overlaps(emvar_regions = TC28_emVars$region, 
                                                active_regions = TC28_active$region, 
                                                filename = "290k_TC28", 
                                                label = "TC28")

K562_overlap_sharing <- MPRA_multi_bed_overlaps(emvar_regions = K562_emVars$region, 
                                                active_regions = K562_active$region, 
                                                filename = "290k_K562", 
                                                label = "K562")

overlap_sharing <-rbind(all_overlap_sharing, CHON_overlap_sharing, K562_overlap_sharing, TC28_overlap_sharing)

#so activity/emvar status is directly proportional to number of tested elements
#plot without cell type
total_active_plot <- ggplot(data = overlap_sharing %>% dplyr::filter(celltype == "ALL"), 
                            aes(x = total_overlaps, y = active_overlaps)) +
  geom_point() + 
  theme_bw() + 
  labs(x = "Total regions per tissue", 
       y = "Regions with activity \nper tissue") + 
  stat_cor(method="pearson", p.accuracy = 1e-7)

unique_active_plot <- ggplot(data = overlap_sharing %>% dplyr::filter(celltype == "ALL"), 
                             aes(x = unique_overlaps, y = active_unique_overlaps)) +
  geom_point() + 
  theme_bw() + 
  labs(x = "Total unique regions per tissue", 
       y = "Unique regions with \nactivity per tissue") + 
  stat_cor(method="pearson", p.accuracy = 1e-7)

#so activity/emvar status is directly proportional to number of tested elements
#plot by cell type
total_active_plot_ct <- ggplot(data = overlap_sharing %>% dplyr::filter(celltype != "ALL"), 
                               aes(x = total_overlaps, y = active_overlaps, shape = celltype)) +
  geom_point() + 
  theme_bw() + 
  labs(x = "Total tested regulatory elements per tissue", 
       y = "Elements with activity per tissue") + 
  stat_cor(method="pearson", p.accuracy = 1e-7)

total_emvar_plot_ct <- ggplot(data = overlap_sharing %>% dplyr::filter(celltype != "ALL"), 
                              aes(x = total_overlaps, y = emvar_overlaps, shape = celltype)) +
  geom_point() + 
  theme_bw() + 
  labs(x = "Total tested regulatory elements per tissue", 
       y = "Elements with differential activity per tissue") + 
  stat_cor(method="pearson", p.accuracy = 1e-7)

unique_active_plot_ct <- ggplot(data = overlap_sharing %>% dplyr::filter(celltype != "ALL"), 
                                aes(x = unique_overlaps, y = active_unique_overlaps, shape = celltype)) +
  geom_point() + 
  theme_bw() + 
  labs(x = "Unique regulatory elements per tissue", 
       y = "Unique elements with activity per tissue") + 
  stat_cor(method="pearson", p.accuracy = 1e-7)

unique_emvar_plot_ct <- ggplot(data = overlap_sharing %>% dplyr::filter(celltype != "ALL"), 
                               aes(x = unique_overlaps, y = emvar_unique_overlaps, shape = celltype))+
  geom_point() + 
  theme_bw() + 
  labs(x = "Unique regulatory elements per tissue", 
       y = "Unique elements with differential activity per tissue") + 
  stat_cor(method="pearson", p.accuracy = 1e-7)


#pleiotropy/timepoint
#overlaps with human signals
human_tp_df_active <- data.frame(celltype = c(rep("Total", 6), rep("CHON002", 6), rep("K562", 6), rep("TC28", 6)), 
                          timepoint = rep(c("Early", "Late", "Shared"), each = 2),
                          class = rep(c("Active", "Inactive"), 6),
                          count = NA)
human_tp_df_active$tp <- factor(human_tp_df_active$timepoint, levels = c("Early", "Late", "Shared"))
human_tp_df_active$celltype <- factor(human_tp_df_active$celltype, levels = c("Total", "CHON002", "K562", "TC28"))

active_regions = unique(c(CHON002_active$region, TC28_active$region, K562_active$region))
active_snps = unique(c(CHON002_active$SNP, TC28_active$SNP, K562_active$SNP))

cartilage_metadata_experimental$activity <- "Inactive"
cartilage_metadata_experimental$activity[which(cartilage_metadata_experimental$region %in% active_regions)] <- "Active"

cartilage_metadata_experimental_regions_simple <- cartilage_metadata_experimental %>% 
  dplyr::select(timepoint, region, activity) %>% unique()

tp_activity_df <- table(cartilage_metadata_experimental_regions_simple$timepoint, cartilage_metadata_experimental_regions_simple$activity) %>% as.data.frame()


#create a function to fill out the dataframe per cell type

#return 6 values, emvars, active only, inactive, for each timepoint
get_pie_tp_data_active <- function(active_df){
  #get data on total library overlaps
  early_regions <- unique(cartilage_metadata_experimental$region[which(cartilage_metadata_experimental$timepoint =="EARLY_ONLY")])
  early_total <- length(early_regions)
  late_regions <- unique(cartilage_metadata_experimental$region[which(cartilage_metadata_experimental$timepoint =="LATE_ONLY")])
  late_total <- length(late_regions)
  shared_regions <- unique(cartilage_metadata_experimental$region[which(cartilage_metadata_experimental$timepoint =="SHARED")])
  shared_total <- length(shared_regions)
  
  #get active 
  early_active <- unique(active_df$region[which(active_df$timepoint == "EARLY_ONLY")])
  late_active <- unique(active_df$region[which(active_df$timepoint == "LATE_ONLY")])
  shared_active <- unique(active_df$region[which(active_df$timepoint == "SHARED")])
  
  return(c(length(early_active), 
           early_total -length(early_active),
           length(late_active),
           late_total -length(late_active),
           length(shared_active),
           shared_total - length(shared_active)))
}

human_tp_df_active$count[1:6] <- get_pie_tp_data_active(rbind(CHON002_active, K562_active, TC28_active))
human_tp_df_active$count[7:12] <- get_pie_tp_data_active(CHON002_active)
human_tp_df_active$count[13:18] <- get_pie_tp_data_active(K562_active)
human_tp_df_active$count[19:24]<- get_pie_tp_data_active(TC28_active)

#plot results
human_tp_celltype_pie_active <- ggplot(data = human_tp_df_active %>% 
                                  group_by(celltype, timepoint) %>% 
                                  mutate(pct = count/sum(count)*100) %>% 
                                  dplyr::filter(celltype != "Total"), aes(x="", y = pct, fill = class)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  #geom_label_repel(aes(label = count, fill = NULL),
  #                 size = 1.5, nudge_x = 1, show.legend = FALSE)+
  geom_text(aes(label = count, x = 1.6,), position = position_stack(vjust=0.5), size = 2) +
  facet_grid(rows = vars(celltype), cols = vars(timepoint), switch = "y") +
  theme_void() +
  theme(legend.position = "bottom") + 
  scale_fill_manual(values = c("#56B4E9", "grey"), name = NULL)

ggsave(plot = human_tp_celltype_pie, filename = "~/Desktop/290k_MPRA_Cluster_Results/results/290k_MPRA_human_timepoint_celltype_activity_pie.png", device = "png", dpi = 300, height = 5, width = 5, units = "in", bg="white")

#count totals for each type
total_tissue_freq <- cartilage_metadata_experimental %>% 
  dplyr::select(tissue_overlap_n) %>% 
  group_by(tissue_overlap_n) %>% 
  summarize(n = n()) %>% mutate(celltype = "ALL")

#get frequency of number of tissues with accessibility per region
all_active_tissue_freq <- cartilage_metadata_experimental %>% 
  dplyr::filter(region %in% unique(c(K562_active$region, TC28_active$region, CHON002_active$region))) %>% 
  dplyr::select(tissue_overlap_n) %>% 
  group_by(tissue_overlap_n) %>% 
  summarize(n = n()) %>% mutate(celltype = "Active")

#generate the plot
total_tissue_v_active_freq_plot <- ggplot(data = total_tissue_freq,
                                          aes(x = tissue_overlap_n, y = n, fill = celltype)) + 
  geom_bar(stat="identity", color = "black") + 
  geom_bar(data = all_active_tissue_freq, stat="identity") + 
  theme_bw() + 
  labs(x = "Number of tissues with accessibility", y  = "Frequency", fill = "") + 
  theme(legend.position = c(0.5, 0.8), 
        legend.direction="horizontal", 
        legend.background = element_rect(fill="lightgray", linewidth = 0.5, linetype="solid")) + 
  scale_fill_manual(breaks=c('ALL', 'Active'), values = c("white", "#56B4E9"), name = NULL)

#now for each cell type

CHON_active_tissue_freq <- CHON002_active %>% 
  dplyr::select(tissue_overlap_n) %>% 
  group_by(tissue_overlap_n) %>% 
  summarize(n = n()) %>% mutate(celltype = "CHON", class = "Active")

TC28_active_tissue_freq <- TC28_active %>% 
  dplyr::select(tissue_overlap_n) %>% 
  group_by(tissue_overlap_n) %>% 
  summarize(n = n()) %>% mutate(celltype = "TC28", class = "Active")

CHON_emvar_tissue_freq <- CHON002_emVars %>% 
  dplyr::select(tissue_overlap_n) %>% 
  group_by(tissue_overlap_n) %>% 
  summarize(n = n()) %>% mutate(celltype = "CHON", class = "emVars")

TC28_emvar_tissue_freq <- TC28_emVars %>% 
  dplyr::select(tissue_overlap_n) %>% 
  group_by(tissue_overlap_n) %>% 
  summarize(n = n()) %>% mutate(celltype = "TC28", class = "emVars")

K562_emvar_tissue_freq <- K562_emVars %>% 
  dplyr::select(tissue_overlap_n) %>% 
  group_by(tissue_overlap_n) %>% 
  summarize(n = n()) %>% mutate(celltype = "K562", class = "emVars")


K562_active_tissue_freq <- K562_active %>% 
  dplyr::select(tissue_overlap_n) %>% 
  group_by(tissue_overlap_n) %>% 
  summarize(n = n()) %>% mutate(celltype = "K562", class = "Active")


tissue_freq_df <- rbind(CHON_active_tissue_freq, 
                        TC28_active_tissue_freq, 
                        K562_active_tissue_freq,
                        CHON_emvar_tissue_freq, 
                        TC28_emvar_tissue_freq, 
                        K562_emvar_tissue_freq)

tissue_freq_df$class[which(tissue_freq_df$class == "emVars")] <- "Differentially Active"
#not sure if these is necessary, has a bug
 tissue_freq_df_norm <- tissue_freq_df %>% 
   pivot_wider(names_from = class, values_from = n) %>% 
   mutate(emvar_perc = `Differentially Active`/Active)
 
 tissue_perc_plot <- ggplot(data = tissue_freq_df_norm, 
                           aes(x = tissue_overlap_n, y = emvar_perc)) + 
   geom_bar(stat="identity") + 
   facet_grid(rows = vars(celltype), scales = "free") +
   theme_bw() + 
   labs(x = "Number of tissues with accessibility", y  = "Proportion of differentially active \nper active sequences")


tissue_freq_plot <- ggplot(data = tissue_freq_df, 
                           aes(x = tissue_overlap_n, y = n)) + 
  geom_bar(stat="identity") + 
  facet_grid(rows = vars(celltype), cols = vars(class), scales = "free") +
  theme_bw() + 
  labs(x = "Number of tissues with accessibility", y  = "Frequency")

tissue_freq_plot_cell_type <- ggdraw() +
  draw_plot(tissue_freq_plot, x = 0, y = .5, width = 1, height = 0.5) +
  draw_plot(tissue_perc_plot, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B"), 
                  size = 12,
                  x = c(0, 0), 
                  y = c(1, 0.5)) + 
  bgcolor("white")

ggsave(plot = tissue_freq_plot_cell_type, filename = "~/Desktop/290k_MPRA_Cluster_Results/results/tissue_freq_plot_cell_type.png", device = "png", dpi = 300, height = 8, width = 6, units = "in", bg="white")

#assemble all the components of a massive active tiles figure

#phyloP for active
active_phylop <- read.table(file = "~/Desktop/290k_MPRA_Cluster_Results/beds/290k_active_tiles_phyloP_hg38.bed", header = F)
active_phylop$Set <- "Active"

inactive_phylop <- read.table(file = "~/Desktop/290k_MPRA_Cluster_Results/beds/290k_inactive_tiles_phyloP_hg38.bed", header = F)
inactive_phylop$Set <- "Inactive"

act_inact_phylop <- rbind(active_phylop, inactive_phylop)
colnames(act_inact_phylop)[5] <- "PhyloP"

wilcox.test(active_phylop$V5, inactive_phylop$V5)

active_phyloP <- ggplot(data = act_inact_phylop, aes(x = Set, y = PhyloP, fill = Set))+
  #geom_violin() +
  geom_boxplot(outliers = F)  + labs(x = NULL, y = "PhyloP") + 
  scale_fill_manual(breaks=c('Inactive', 'Active'), values = c("grey", "#56B4E9")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  annotate("text", label='NS', x = 1.5, y = 2.5) + 
  ylim(0, 2.5)

#create beds
cartilage_metadata_experimental %>% 
  dplyr::filter(region %in% unique(c(CHON002_active$region, TC28_active$region, K562_active$region))) %>% 
  mutate(region_end = region_start + 540) %>% 
  dplyr::select(chr, region_start, region_end) %>% 
  unique() %>% 
  write_bed(filename = "~/Desktop/290k_MPRA_Cluster_Results/beds/290k_active_regions_hg38.bed")

cartilage_metadata_experimental %>% 
  dplyr::filter(!region %in% unique(c(CHON002_active$region, TC28_active$region, K562_active$region))) %>% 
  mutate(region_end = region_start + 540) %>% 
  dplyr::select(chr, region_start, region_end) %>% 
  unique() %>% 
  write_bed(filename = "~/Desktop/290k_MPRA_Cluster_Results/beds/290k_inactive_regions_hg38.bed")

#TSS for active
active_TSS <- read.table(file = "~/Desktop/290k_MPRA_Cluster_Results/beds/290k_active_regions_TSS_hg38.bed", header = F)
active_TSS$Set <- "Active"

inactive_TSS <- read.table(file = "~/Desktop/290k_MPRA_Cluster_Results/beds/290k_inactive_regions_TSS_hg38.bed", header = F)
inactive_TSS$Set <- "Inactive"

act_inact_TSS <- rbind(active_TSS, inactive_TSS)
colnames(act_inact_TSS)[10] <- "TSS"
act_inact_TSS$kTSS <- act_inact_TSS$TSS/1000

wilcox.test(active_TSS$V10, inactive_TSS$V10)

active_TSS_plot <- ggplot(data = act_inact_TSS, aes(x = Set, y = kTSS, fill = Set))+
  #geom_violin() +
  geom_boxplot(outliers = F)  + 
  labs(x = NULL, y = "Distance to TSS\n(Kilobases)") + 
  scale_fill_manual(breaks=c('Inactive', 'Active'), values = c("grey", "#56B4E9")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  annotate("text", label='***', x = 1.5, y = 45) + 
  ylim(0, 45)

human_tp_total_pie_active <- ggplot(data = human_tp_df_active %>% 
                               group_by(celltype, timepoint) %>% 
                               mutate(pct = count/sum(count)*100) %>% 
                               dplyr::filter(celltype == "Total"), aes(x="", y = pct, fill = class)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  #geom_label_repel(aes(label = count, fill = NULL),
  #                 size = 1.5, nudge_x = 1, show.legend = FALSE)+
  geom_text(aes(label = formatC(count, format = "d", big.mark = ","), x = 1), position = position_stack(vjust=0.5), size = 3) +
  facet_grid(cols = vars(timepoint), switch = "x") +
  theme_void() +
  theme(legend.position = "right", 
        strip.placement = "outside", 
        plot.margin = unit(c(0, 0, 0.5, 0), "cm")) + 
  scale_fill_manual(values = c("#56B4E9", "grey"), name = NULL) 

active_total_plot <- ggdraw() +
  draw_plot(activity_pie, x = 0, y = .7, width = 0.4, height = .3) +
  draw_plot(activity_venn, x = 0.4, y = .7, width = 0.6, height = .3) +
  draw_plot(total_active_plot, x = 0, y = .45, width = 0.5, height = .25) +
  draw_plot(unique_active_plot, x = 0.5, y = .45, width = 0.5, height = .25) +
  draw_plot(total_tissue_v_active_freq_plot, x = 0, y = 0.2, width = 0.5, height = .25) +
  draw_plot(active_phyloP, x = 0.5, y = 0.2, width = 0.25, height = .25) +
  draw_plot(active_TSS_plot, x = 0.75, y = 0.2, width = 0.25, height = .25) +
  draw_plot(human_tp_total_pie_active, x = 0, y = 0, width = 1, height = .2) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F", "G", "H"), 
                  size = 10,
                  x = c(0, 0.4, 0, 0.5, 0, 0.5, 0.75, 0), 
                  y = c(1, 1, 0.7, 0.7, 0.45, 0.45, 0.45, 0.2)) + 
  bgcolor("white")

ggsave(plot = active_total_plot + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Cartilage MPRA Paper/Code/fig2_activity_plot.png", 
       device = "png", dpi = 300, height = 9.5, width = 6.5, units = "in")

#motif enrichment in different sets

#load library to try TF analysis
library("monaLisa")
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(mia)

#get pwm matrix for TFs
pwms <- getMatrixSet(JASPAR2020,
                     opts = list(matrixtype = "PWM",
                                 species    = "9606"))

pairwise_TF_analysis_MPRA <- function(sequences1, 
                                      sequences2, 
                                      ids =c("bed1", "bed2")){
  seqs1 <- DNAStringSet(sequences1)
  seqs2 <- DNAStringSet(sequences2)
  
  #combine sequences into single object
  comp_seqs <-  c(seqs1, seqs2)
  bins <- rep(ids, c(length(seqs1), length(seqs2)))
  bins <- factor(bins)
  #perform motif enrichment
  se <- calcBinnedMotifEnrR(seqs = comp_seqs, bins = bins,
                            pwmL = pwms, BPPARAM = BiocParallel::MulticoreParam(1))
  se_long <- merge(meltAssay(se, assay.type = "negLog10Padj", add_row_data = T), meltAssay(se, assay.type = "log2enr", add_row_data = T))
  
  return(se_long)
}

CHON_TFs <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% CHON002_active$SNP)], 
                                      sequences2 = cartilage_metadata_experimental$human, 
                                      ids = c("CHON_active", "all") )%>% 
  as_tibble() %>% dplyr::filter(negLog10Padj > 4)


TC28_TFs <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% TC28_active$SNP)], 
                                      sequences2 = cartilage_metadata_experimental$human, 
                                      ids = c("TC28_active", "all")
) %>% as_tibble() %>% dplyr::filter(negLog10Padj > 4)

K562_TFs <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% K562_active$SNP)], 
                                      sequences2 = cartilage_metadata_experimental$human, 
                                      ids = c("K562_active", "all")
) %>% as_tibble() %>% dplyr::filter(negLog10Padj > 4)

depleted_TFs <- intersect(CHON_TFs$motif.name[which(CHON_TFs$SampleID == "all")], intersect(TC28_TFs$motif.name[which(TC28_TFs$SampleID == "all")], K562_TFs$motif.name[which(K562_TFs$SampleID == "all")]))

# make Venn Diagrams
activity_human_TF_venn <- ggVennDiagram(
  x = list(unique(CHON_TFs$motif.name), unique(K562_TFs$motif.name), unique(TC28_TFs$motif.name)),
  category.names = c("  CHON002" , "K562" , "TC28") 
) + labs(fill = "TF Count")

#get unique TFs
chon_active_unique_TFs <- CHON_TFs %>%  dplyr::filter(motif.name %in% setdiff(CHON_TFs$motif.name, c(TC28_TFs$motif.name, K562_TFs$motif.name)))
tc28_active_unique_TFs <- TC28_TFs %>%  dplyr::filter(motif.name %in% setdiff(TC28_TFs$motif.name, c(CHON_TFs$motif.name, K562_TFs$motif.name)))
k562_active_unique_TFs <- K562_TFs %>%  dplyr::filter(motif.name %in% setdiff(K562_TFs$motif.name, c(TC28_TFs$motif.name, CHON_TFs$motif.name)))


#repeat for chimp
CHON_TFs_chimp <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% CHON002_active$SNP)], 
                                            sequences2 = cartilage_metadata_experimental$chimp, 
                                            ids = c("CHON_active", "all")
) %>% as_tibble() %>% dplyr::filter(log2enr > 0 & negLog10Padj > 4 & SampleID == "CHON_active")


TC28_TFs_chimp <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% TC28_active$SNP)], 
                                            sequences2 = cartilage_metadata_experimental$chimp, 
                                            ids = c("TC28_active", "all")
) %>% as_tibble() %>% dplyr::filter(log2enr > 0 & negLog10Padj > 4 & SampleID == "TC28_active")

K562_TFs_chimp <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% K562_active$SNP)], 
                                            sequences2 = cartilage_metadata_experimental$chimp, 
                                            ids = c("K562_active", "all")
) %>% as_tibble() %>% dplyr::filter(log2enr > 0 & negLog10Padj > 4 & SampleID == "K562_active")

# make Venn Diagrams
activity_chimp_TF_venn <- ggVennDiagram(
  x = list(unique(CHON_TFs_chimp$motif.name), unique(K562_TFs_chimp$motif.name), unique(TC28_TFs_chimp$motif.name)),
  category.names = c("  CHON002" , "K562" , "TC28") 
) + labs(fill = "TF Count")

#get unique TFs
chon_active_unique_TFs_chimp <- CHON_TFs_chimp %>%  dplyr::filter(motif.name %in% setdiff(CHON_TFs_chimp$motif.name, c(TC28_TFs_chimp$motif.name, K562_TFs_chimp$motif.name)))
tc28_active_unique_TFs_chimp <- TC28_TFs_chimp %>%  dplyr::filter(motif.name %in% setdiff(TC28_TFs_chimp$motif.name, c(CHON_TFs_chimp$motif.name, K562_TFs_chimp$motif.name)))
k562_active_unique_TFs_chimp <- K562_TFs_chimp %>%  dplyr::filter(motif.name %in% setdiff(K562_TFs_chimp$motif.name, c(TC28_TFs_chimp$motif.name, CHON_TFs_chimp$motif.name)))

#create active TF tables for supplement
TC28_TFs$celltype <- NULL
CHON_TFs$celltype <- NULL
K562_TFs$celltype <- NULL

human_TFs <- rbind(CHON_TFs, TC28_TFs, K562_TFs) %>% dplyr::select(-motif.pfm, -motif.pwm)
write_delim(human_TFs, delim = "\t", file = "/Users/alexanderokamoto/Desktop/290k_MPRA_Cluster_Results/results/human_active_TFs.txt")

chimp_TFs <- rbind(CHON_TFs_chimp, TC28_TFs_chimp, K562_TFs_chimp) %>% dplyr::select(-motif.pfm, -motif.pwm)
write_delim(chimp_TFs, delim = "\t", file = "/Users/alexanderokamoto/Desktop/290k_MPRA_Cluster_Results/results/chimp_active_TFs.txt")

active_TF_plot <- ggdraw() +
  draw_plot(activity_human_TF_venn, x = 0, y = 0.5, width = 1, height = .5) +
  draw_plot(activity_chimp_TF_venn, x = 0, y = 0, width = 1, height = .5) +
  draw_plot_label(label = c("Human", "Chimp"), 
                  size = 10,
                  x = c(0, 0), 
                  y = c(1, 0.5)) + 
  bgcolor("white")

ggsave(plot = active_TF_plot + panel_border(color = "black", size = 1), 
       filename = "active_TF_enrichment_plot.png", 
       device = "png", dpi = 300, height = 6., width = 3, units = "in")

#now direct species comparison
#repeat for chimp
CHON_TFs_human_chimp <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% CHON002_active$SNP)], 
                                                  sequences2 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% CHON002_active$SNP)], 
                                                  ids = c("CHON_active_human", "CHON_active_chimp")
) %>% as_tibble() %>% dplyr::filter(log2enr > 0 & negLog10Padj > 4)


TC28_TFs_human_chimp <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% TC28_active$SNP)], 
                                                  sequences2 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% TC28_active$SNP)], 
                                                  ids = c("TC28_active_human", "TC28_active_chimp")
) %>% as_tibble() %>% dplyr::filter(log2enr > 0 & negLog10Padj > 4)

K562_TFs_human_chimp <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% K562_active$SNP)], 
                                                  sequences2 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% K562_active$SNP)], 
                                                  ids = c("K562_active_human", "K562_active_chimp")
) %>% as_tibble() %>% dplyr::filter(log2enr > 0 & negLog10Padj > 4)

# make Venn Diagrams
activity_chimp_TF_venn <- ggVennDiagram(
  x = list(unique(CHON_TFs_chimp$motif.name), unique(K562_TFs_chimp$motif.name), unique(TC28_TFs_chimp$motif.name)),
  category.names = c("  CHON002" , "K562" , "TC28") 
) + 
  labs(fill = "TF Count")


#add new activity figure 


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

active_plotting_df <- emvars_df %>% 
  dplyr::select(celltype, active_controls, neg_controls) %>% pivot_longer(!celltype, names_to = "Set", values_to = "percent")
active_plotting_df <- rbind(active_plotting_df, c("CHON002", "Experimental", CHON_active_regions_count/total_regions_counts*100))
active_plotting_df <- rbind(active_plotting_df, c("TC28", "Experimental", TC28_active_regions_count/total_regions_counts*100))
active_plotting_df <- rbind(active_plotting_df, c("K562", "Experimental", K562_active_regions_count/total_regions_counts*100))
active_plotting_df$percent <- as.numeric(active_plotting_df$percent)

#clean up set labels for plotting
active_plotting_df$Set <- gsub(pattern = "active_controls", x = active_plotting_df$Set, replacement = "Positive Control")
active_plotting_df$Set <- gsub(pattern = "neg_controls", x = active_plotting_df$Set, replacement = "Negative Control")

active_plotting_df$Set <- factor(active_plotting_df$Set, levels = c("Negative Control", "Positive Control", "Experimental"))

activity_plot <- ggplot(active_plotting_df, aes(x = Set, y = percent, fill = Set)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~celltype) + 
  labs(x = NULL, y = "Percent Active", fill = NULL) + 
  theme_classic() + 
  scale_fill_manual(values = c("#CECECE", "#7D7D7D","#000000"))+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "bottom")

active_plotting_df <- emvars_df %>% 
  dplyr::select(celltype, active_controls, neg_controls) %>% pivot_longer(!celltype, names_to = "Set", values_to = "percent")
active_plotting_df <- rbind(active_plotting_df, c("CHON002", "Experimental", CHON_active_regions_count/total_regions_counts*100))
active_plotting_df <- rbind(active_plotting_df, c("TC28", "Experimental", TC28_active_regions_count/total_regions_counts*100))
active_plotting_df <- rbind(active_plotting_df, c("K562", "Experimental", K562_active_regions_count/total_regions_counts*100))
active_plotting_df$percent <- as.numeric(active_plotting_df$percent)


emvar_plotting_df <- data.frame(
  celltype = c("CHON002", "K562", "TC28"), 
  Percent = c(length(unique(CHON002_emVars$region))/CHON_active_regions_count*100, 
              length(unique(K562_emVars$region))/K562_active_regions_count*100,
              length(unique(TC28_emVars$region))/TC28_active_regions_count*100)
  )

emvar_plot <- ggplot(emvar_plotting_df, aes(x = celltype, y = Percent)) + 
  geom_bar(stat = "identity") + 
  labs(x = NULL, y = "Percent Differentally Active") + 
  theme_classic() # + 
  #scale_fill_manual(values = c("#CECECE", "#7D7D7D","#000000"))+
  #theme(axis.text.x = element_blank(), 
  #      axis.ticks.x = element_blank())

activity_summary <- ggdraw() +
  draw_plot(activity_plot , x = 0, y = 0, width = 0.65, height = 1) +
  draw_plot(emvar_plot, x = 0.65, y = 0, width = 0.35, height = 1) +
  draw_plot_label(label = c("A", "B"), 
                  size = 10,
                  x = c(0, 0.65), 
                  y = c(1, 1))  + 
  bgcolor("white")

ggsave(plot = activity_summary , 
       filename = "~/Dropbox/Cartilage MPRA Paper/Code/fig2_activity_summary_plot.png", 
       device = "png", dpi = 300, height = 3.5, width = 6.5, units = "in")

