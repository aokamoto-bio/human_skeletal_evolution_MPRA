## Code for Massively parallel functional screen identifies thousands of regulatory differences in human versus chimpanzee postcranial skeletal development

## Part 4: Investigation of evolutionary signals in MPRA results

#load data
#modify to specify file path for your computer
setwd("~/Desktop/290k_MPRA_Cluster_Results/results/")
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggVennDiagram)
library(cowplot)
library(magick) # for adding images to ggplots
library(ggbreak)

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

#overlaps with human signals
human_signal_df <- data.frame(celltype = c(rep("Total", 9), rep("CHON002", 9), rep("K562", 9), rep("TC28", 9)), 
                              signal = rep(c("HARs", "HAQERs", "hCONDELs"), each = 3),
                              class = rep(c("Differentially\nActive", "Active", "Inactive"), 12),
                              count = NA)
human_signal_df$signal <- factor(human_signal_df$signal, levels = c("HARs", "HAQERs", "hCONDELs"))
human_signal_df$celltype <- factor(human_signal_df$celltype, levels = c("Total", "CHON002", "K562", "TC28"))

#create a function to fill out the dataframe per cell type
#return 9 values, emvars, active only, inactive, for HAR, HAQER, hCONDEL
get_pie_data <- function(active_df, emvars_df){
  #get data on total library overlaps
  HAR_regions <- unique(cartilage_metadata_experimental$region[which(!is.na(cartilage_metadata_experimental$HARs_Overlapped))])
  HAR_total <- length(HAR_regions)
  HAQER_regions <- unique(cartilage_metadata_experimental$region[which(!is.na(cartilage_metadata_experimental$HAQERs_Overlapped))])
  HAQER_total <- length(HAQER_regions)
  hCONDEL_regions <- unique(cartilage_metadata_experimental$region[which(!is.na(cartilage_metadata_experimental$hCONDELs_Overlapped))])
  hCONDEL_total <- length(hCONDEL_regions)
  
  #get active 
  HAR_active <- unique(active_df$region[which(!is.na(active_df$HARs_Overlapped))])
  HAQER_active <- unique(active_df$region[which(!is.na(active_df$HAQERs_Overlapped))])
  hCONDEL_active <- unique(active_df$region[which(!is.na(active_df$hCONDELs_Overlapped))])
  
  #get active 
  HAR_emVars <- unique(emvars_df$region[which(!is.na(emvars_df$HARs_Overlapped))])
  HAQER_emVars <- unique(emvars_df$region[which(!is.na(emvars_df$HAQERs_Overlapped))])
  hCONDEL_emVars <- unique(emvars_df$region[which(!is.na(emvars_df$hCONDELs_Overlapped))])
  
  return(c(length(HAR_emVars), 
           length(HAR_active) - length(HAR_emVars),
           HAR_total - length(HAR_active),
           length(HAQER_emVars), 
           length(HAQER_active) - length(HAQER_emVars),
           HAQER_total -length(HAQER_active),
           length(hCONDEL_emVars), 
           length(hCONDEL_active) - length(hCONDEL_emVars),
           hCONDEL_total -length(hCONDEL_active)))
}

human_signal_df$count[1:9] <- get_pie_data(rbind(CHON002_active, K562_active, TC28_active), 
                                           emvars_df = rbind(CHON002_emVars, K562_emVars, TC28_emVars))
human_signal_df$count[10:18] <- get_pie_data(CHON002_active, emvars_df = CHON002_emVars)
human_signal_df$count[19:27] <- get_pie_data(K562_active, emvars_df = K562_emVars)
human_signal_df$count[28:36]<- get_pie_data(TC28_active, emvars_df = TC28_emVars)
human_signal_df <- human_signal_df %>% dplyr::filter(signal != "hCONDELs")

human_signal_celltype_pie <- ggplot(data = human_signal_df %>% 
                                      group_by(celltype, signal) %>% 
                                      mutate(pct = count/sum(count)*100), aes(x="", y = pct, fill = class)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label = count, x = 1.1,), position = position_stack(vjust=0.5), size = 3) +
  facet_grid(cols = vars(signal), switch = "y") +
  theme_void() +
  theme(legend.position = "bottom") + 
  scale_fill_manual(values = c("#56B4E9", "#0079B2", "grey"), name = NULL)

ggsave(plot = human_signal_celltype_pie, filename = "290k_MPRA_human_signal_celltype_pie.png", device = "png", dpi = 300, height = 7, width = 5, units = "in", bg="white")

human_signal_df$class <- factor(human_signal_df$class, levels = c("Differentially\nActive", "Active", "Inactive"))
#smaller plot for figure 
HAR_HAQER_pie <- ggplot(data = human_signal_df %>% 
            group_by(celltype, signal) %>% 
            dplyr::filter(celltype == "Total") %>% 
            dplyr::filter(class != "Inactive") %>% 
            mutate(pct = count/sum(count)*100), aes(x="", y = pct, fill = class)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label = count, x = 1.1,), position = position_stack(vjust=0.5), size = 3) +
  facet_grid(cols = vars(signal), switch = "y") +
  theme_void() +
  theme(legend.position = "bottom") + 
  scale_fill_manual(values = c("#0079B2", "#56B4E9"), 
                    name = NULL, labels = c("Differentially Active", "Active")) +
  labs(x =NULL)

#emvar versus hars/haqer/hcondel 4x3 (regular, unique)
emvar_dfs <- list(CHON002_emVars, TC28_emVars, K562_emVars)
CHON_emvar_bed <- get_MPRA_bed(unique(CHON002_emVars$region), filename = "290k_CHON002_emvars")
CHON_emvar_overlaps <- bed_capellini_skeleton_overlaps(bed_file = CHON_emvar_bed, label = "CHON_emVars")
TC28_emvar_bed <- get_MPRA_bed(TC28_emVars$region, filename = "290k_TC28_emvars")
bed_capellini_skeleton_overlaps(bed_file = TC28_emvar_bed, label = "TC28_emVars")
K562_emvar_bed <- get_MPRA_bed(K562_emVars$region, filename = "290k_K562_emvars")
bed_capellini_skeleton_overlaps(bed_file = K562_emvar_bed, label = "K562_emVars")
all_emvar_bed <- get_MPRA_bed(unique(c(CHON002_emVars$region, TC28_emVars$region, K562_emVars$region)), filename = "290k_all_emvars")
cell_emvar_beds <- c(all_emvar_bed, CHON_emvar_bed, TC28_emvar_bed, K562_emvar_bed)
cell_emvar_beds_HARs <- gsub(x = cell_emvar_beds, pattern = "hg38.bed", replacement = "HARs_hg38.bed")
cell_emvar_beds_HAQERs <- gsub(x = cell_emvar_beds, pattern = "hg38.bed", replacement = "HAQERs_hg38.bed")

cell_names  <- c("all", "CHON", "TC28", "K562")

for(c in 1:3){
  #hars
  emvar_dfs[c] %>% as.data.frame() %>% 
    dplyr::filter(!is.na(HARs_Overlapped)) %>% 
    mutate(region_end = (region_start+540)) %>% 
    dplyr::select(chr.x, region_start, region_end, HARs_Overlapped) %>% 
    write_bed(file = cell_emvar_beds_HARs[c+1], ncol = 4)
  
  #HAQERS
  emvar_dfs[c] %>% as.data.frame() %>% 
    dplyr::filter(!is.na(HAQERs_Overlapped)) %>% 
    mutate(region_end = (region_start+540)) %>% 
    dplyr::select(chr.x, region_start, region_end, HAQERs_Overlapped) %>% 
    write_bed(file = cell_emvar_beds_HAQERs[c+1], ncol = 4)
}

system(paste("cat", paste(cell_emvar_beds_HARs[2:4], collapse = " "), "| sort -k 1,1 -k2,2n | bedtools merge -i stdin >", cell_emvar_beds_HARs[1], sep = " "))
system(paste("cat", paste(cell_emvar_beds_HAQERs[2:4], collapse = " "), "| sort -k 1,1 -k2,2n | bedtools merge -i stdin >", cell_emvar_beds_HAQERs[1], sep = " "))

for(c in 1:4){
  bed_capellini_skeleton_overlaps(bed_file = cell_emvar_beds_HARs[c], 
                                  label = paste(cell_names[c], "_emvar_HAR", sep = ""))
  bed_capellini_skeleton_overlaps(bed_file = cell_emvar_beds_HAQERs[c], 
                                  label = paste(cell_names[c], "_emvar_HAQER", sep = ""))
  bed_capellini_skeleton_overlaps(bed_file = cell_emvar_beds_HARs[c], 
                                  label = paste(cell_names[c], "_emvar_HAR", sep = ""), 
                                  timepoint = "unique")
  bed_capellini_skeleton_overlaps(bed_file = cell_emvar_beds_HAQERs[c], 
                                  label = paste(cell_names[c], "_emvar_HAQER", sep = ""), 
                                  timepoint = "unique")
} #end for loop

#plot results

all_emvar_HARs_skeleton <- image_read("~/Desktop/Autopod_Skeleton_Figures/all_emvar_HAR_both_skeleton_heatmap.png") %>%
  image_ggplot()
all_emvar_HAQERs_skeleton <- image_read("~/Desktop/Autopod_Skeleton_Figures/all_emvar_HAQER_both_skeleton_heatmap.png") %>%
  image_ggplot()

CHON_emvar_HARs_skeleton <- image_read("~/Desktop/Autopod_Skeleton_Figures/CHON_emvar_HAR_both_skeleton_heatmap.png") %>%
  image_ggplot()
CHON_emvar_HAQERs_skeleton <- image_read("~/Desktop/Autopod_Skeleton_Figures/CHON_emvar_HAQER_both_skeleton_heatmap.png") %>%
  image_ggplot()

TC28_emvar_HARs_skeleton <- image_read("~/Desktop/Autopod_Skeleton_Figures/TC28_emvar_HAR_both_skeleton_heatmap.png") %>%
  image_ggplot()
TC28_emvar_HAQERs_skeleton <- image_read("~/Desktop/Autopod_Skeleton_Figures/TC28_emvar_HAQER_both_skeleton_heatmap.png") %>%
  image_ggplot()

K562_emvar_HARs_skeleton <- image_read("~/Desktop/Autopod_Skeleton_Figures/K562_emvar_HAR_both_skeleton_heatmap.png") %>%
  image_ggplot()
K562_emvar_HAQERs_skeleton <- image_read("~/Desktop/Autopod_Skeleton_Figures/K562_emvar_HAQER_both_skeleton_heatmap.png") %>%
  image_ggplot()

#make plot by cell type
emvar_HARplus_celltype_overlaps_skeleton <- ggdraw() +
  draw_plot(all_emvar_HARs_skeleton, x = 0, y = 0.75, width = 0.5, height = 0.25) +
  draw_plot(all_emvar_HAQERs_skeleton, x = 0.5, y = 0.75, width = 0.5, height = 0.25) +
  draw_plot(CHON_emvar_HARs_skeleton, x = 0, y = 0.5, width = 0.5, height = 0.25) +
  draw_plot(CHON_emvar_HAQERs_skeleton, x = 0.5, y = 0.5, width = 0.5, height = 0.25) +
  draw_plot(TC28_emvar_HARs_skeleton, x = 0, y = 0.25, width = 0.5, height = 0.25) +
  draw_plot(TC28_emvar_HAQERs_skeleton, x = 0.5, y = 0.25, width = 0.5, height = 0.25) +
  draw_plot(K562_emvar_HARs_skeleton, x = 0, y = 0, width = 0.5, height = 0.25) +
  draw_plot(K562_emvar_HAQERs_skeleton, x = 0.5, y = 0, width = 0.5, height = 0.25) +
  draw_plot_label(label = c("A", "B", "C", "D"), 
                  size = 12,
                  x = c(0, 0, 0, 0), 
                  y = c(1, 0.75, 0.5, 0.25)) + 
  bgcolor("white")

ggsave(plot = emvar_HARplus_celltype_overlaps_skeleton + panel_border(color = "black", size = 1), 
       filename = "~/Desktop/290k_MPRA_Cluster_Results/results/emvar_HARplus_celltype_overlaps_skeleton.png", 
       device = "png", dpi = 300, height = 9, width = 9, units = "in")

#plot unique results

all_emvar_HARs_skeleton_unique <- image_read("~/Desktop/Autopod_Skeleton_Figures/all_emvar_HAR_unique_skeleton_heatmap.png") %>%
  image_ggplot()
all_emvar_HAQERs_skeleton_unique <- image_read("~/Desktop/Autopod_Skeleton_Figures/all_emvar_HAQER_unique_skeleton_heatmap.png") %>%
  image_ggplot()

CHON_emvar_HARs_skeleton_unique <- image_read("~/Desktop/Autopod_Skeleton_Figures/CHON_emvar_HAR_unique_skeleton_heatmap.png") %>%
  image_ggplot()
CHON_emvar_HAQERs_skeleton_unique <- image_read("~/Desktop/Autopod_Skeleton_Figures/CHON_emvar_HAQER_unique_skeleton_heatmap.png") %>%
  image_ggplot()

TC28_emvar_HARs_skeleton_unique <- image_read("~/Desktop/Autopod_Skeleton_Figures/TC28_emvar_HAR_unique_skeleton_heatmap.png") %>%
  image_ggplot()
TC28_emvar_HAQERs_skeleton_unique <- image_read("~/Desktop/Autopod_Skeleton_Figures/TC28_emvar_HAQER_unique_skeleton_heatmap.png") %>%
  image_ggplot()

K562_emvar_HARs_skeleton_unique <- image_read("~/Desktop/Autopod_Skeleton_Figures/K562_emvar_HAR_unique_skeleton_heatmap.png") %>%
  image_ggplot()
K562_emvar_HAQERs_skeleton_unique <- image_read("~/Desktop/Autopod_Skeleton_Figures/K562_emvar_HAQER_unique_skeleton_heatmap.png") %>%
  image_ggplot()

#make plot by cell type
emvar_HARplus_celltype_overlaps_skeleton_unique <- ggdraw() +
  draw_plot(all_emvar_HARs_skeleton_unique, x = 0, y = 0.75, width = 0.5, height = 0.25) +
  draw_plot(all_emvar_HAQERs_skeleton_unique, x = 0.5, y = 0.75, width = 0.5, height = 0.25) +
  draw_plot(CHON_emvar_HARs_skeleton_unique, x = 0, y = 0.5, width = 0.5, height = 0.25) +
  draw_plot(CHON_emvar_HAQERs_skeleton_unique, x = 0.5, y = 0.5, width = 0.5, height = 0.25) +
  draw_plot(TC28_emvar_HARs_skeleton_unique, x = 0, y = 0.25, width = 0.5, height = 0.25) +
  draw_plot(TC28_emvar_HAQERs_skeleton_unique, x = 0.5, y = 0.25, width = 0.5, height = 0.25) +
  draw_plot(K562_emvar_HARs_skeleton_unique, x = 0, y = 0, width = 0.5, height = 0.25) +
  draw_plot(K562_emvar_HAQERs_skeleton_unique, x = 0.5, y = 0, width = 0.5, height = 0.25) +
  draw_plot_label(label = c("A", "B", "C", "D"), 
                  size = 12,
                  x = c(0, 0, 0, 0), 
                  y = c(1, 0.75, 0.5, 0.25)) + 
  bgcolor("white")

ggsave(plot = emvar_HARplus_celltype_overlaps_skeleton_unique + panel_border(color = "black", size = 1), 
       filename = "~/Desktop/290k_MPRA_Cluster_Results/results/emvar_HARplus_celltype_overlaps_skeleton_unique.png", 
       device = "png", dpi = 300, height = 9, width = 9, units = "in")

#number of bps differences versus effect size

#create dataframes of active and emvar sequences only
all_emvars_df <- cartilage_metadata_experimental %>% 
  dplyr::filter(SNP %in% unique(c(CHON002_emVars$SNP, TC28_emVars$SNP, K562_emVars$SNP)))
all_active_df <- cartilage_metadata_experimental %>% 
  dplyr::filter(SNP %in% unique(c(CHON002_active$SNP, TC28_active$SNP, K562_active$SNP)))

all_emvars_df$CHON_skew <- NA
all_emvars_df$CHON_skew[which(all_emvars_df$SNP %in% CHON002_emVars$SNP)] <- CHON002_emVars$Log2Skew
all_emvars_df$K562_skew <- NA
all_emvars_df$K562_skew[which(all_emvars_df$SNP %in% K562_emVars$SNP)] <- K562_emVars$Log2Skew
all_emvars_df$TC28_skew <- NA
all_emvars_df$TC28_skew[which(all_emvars_df$SNP %in% TC28_emVars$SNP)] <- TC28_emVars$Log2Skew

all_emvars_df$skew_sign_match <- TRUE
all_emvars_df %>% rowwise() %>% mutate(mismatch = length(unique(na.omit(c(sign(CHON_skew), sign(TC28_skew), sign(K562_skew)))))) %>% pull(mismatch) %>% table()
sign_mismatch <- all_emvars_df %>% rowwise() %>% mutate(mismatch = length(unique(na.omit(c(sign(CHON_skew), sign(TC28_skew), sign(K562_skew)))))) %>% dplyr::filter(mismatch == 2) %>% pull(SNP)
all_emvars_df$skew_sign_match[which(all_emvars_df$SNP %in% sign_mismatch)] <- FALSE

#calculate maximum skew in any cell type if consistent sign
all_emvars_df$max_skew <- NA
all_emvars_df$max_skew[which(all_emvars_df$skew_sign_match)] <- all_emvars_df %>% 
  dplyr::filter(skew_sign_match) %>% 
  rowwise() %>% 
  mutate(max_skew2 = max(na.omit(c(abs(CHON_skew), abs(TC28_skew), abs(K562_skew))))) %>% pull(max_skew2)

#add information about HARs, HAQERs, etc.
all_emvars_df$status <- "none"
all_emvars_df$status[which(!is.na(all_emvars_df$HARs_Overlapped))] <- "HAR"
all_emvars_df$status[which(!is.na(all_emvars_df$HAQERs_Overlapped))] <- "HAQER"
all_emvars_df$status[which(!is.na(all_emvars_df$hCONDELs_Overlapped))] <- "hCONDEL"
all_emvars_df$status[which(!is.na(all_emvars_df$hCONDELs_Overlapped)&!is.na(all_emvars_df$HARs_Overlapped))] <- "HAR/hCONDEL"
all_emvars_df$status[which(!is.na(all_emvars_df$HAQERs_Overlapped)&!is.na(all_emvars_df$HARs_Overlapped))] <- "HAR/HAQER"
all_emvars_df$status[which(!is.na(all_emvars_df$HAQERs_Overlapped)&!is.na(all_emvars_df$hCONDELs_Overlapped))] <- "hCONDEL/HAQER"
all_emvars_df$status[which(!is.na(all_emvars_df$HAQERs_Overlapped)&!is.na(all_emvars_df$HARs_Overlapped)&!is.na(all_emvars_df$hCONDELs_Overlapped))] <- "All_Three"

#calculate enrichment of HARs/emvars enrichment
all_active_df_reduced <- all_active_df %>% dplyr::select(region, HARs_Overlapped, HAQERs_Overlapped) %>% unique()
all_active_df_reduced$emvar <- FALSE
all_active_df_reduced$emvar[which(all_active_df_reduced$region %in% all_emvars_df$region)] <- TRUE

emvar_HARs_table <- table(!is.na(all_active_df_reduced$HARs_Overlapped), all_active_df_reduced$emvar)
fisher.test(emvar_HARs_table)

emvar_HAQERs_table <- table(!is.na(all_active_df_reduced$HAQERs_Overlapped), all_active_df_reduced$emvar)
fisher.test(emvar_HAQERs_table)

#first check if within active HAR/HAR controls, there is more diff. activity in HARs
HAR_active_df <- all_active_df %>% dplyr::filter(!is.na(HARs_Overlapped))%>% 
  mutate(type = "HAR") %>% dplyr::select(region, HARs_Overlapped) %>% unique()

HAR_active_controls_df <- all_active_df %>% 
  dplyr::filter(is.na(HARs_Overlapped)) %>% 
  dplyr::filter((bp_diffs > 4 & bp_diffs < 8) |
                  (align_type == "chimp_deletion" & deletion_size < 5)) %>% 
  mutate(type = "CONTROL") %>% dplyr::select(region, HARs_Overlapped) %>% unique()

HAR_active_comp_df <- rbind(HAR_active_df, HAR_active_controls_df)
#add emvar information
HAR_active_comp_df$emvar <- "Active"
HAR_active_comp_df$emvar[which(HAR_active_comp_df $region %in% all_emvars_df$region)] <- "Diff. Active"

HAR_active_comp_table <- table(!is.na(HAR_active_comp_df$HARs_Overlapped), HAR_active_comp_df$emvar)
fisher.test(HAR_active_comp_table)

#now to test for max skew difference
HAR_emvars_df <- all_emvars_df %>% 
  dplyr::filter(!is.na(HARs_Overlapped))%>% 
  mutate(type = "HAR") %>% 
  dplyr::select(region, max_skew) %>% 
  group_by(region) %>% 
  mutate(overall_max_skew = max(max_skew)) %>% 
  ungroup %>% 
  dplyr::select(region, overall_max_skew) %>%
  unique()
HAR_emvars_df$Overlap = "HAR"
HAR_emvars_df$comp = "HAR"

HAR_controls_df <- all_emvars_df %>% 
  dplyr::filter(is.na(HARs_Overlapped)) %>% 
  dplyr::filter((bp_diffs > 4 & bp_diffs < 8) |
                  (align_type == "chimp_deletion" & deletion_size < 5)) %>% 
  mutate(type = "CONTROL") %>% dplyr::select(region, max_skew) %>% 
  group_by(region) %>% 
  mutate(overall_max_skew = max(max_skew)) %>% 
  ungroup %>% 
  dplyr::select(region, overall_max_skew) %>%
  unique()
HAR_controls_df$Overlap = "CTRL"
HAR_controls_df$comp = "HAR"

HAR_comp_df <- rbind(HAR_emvars_df, HAR_controls_df)

wilcox.test(HAR_emvars_df$overall_max_skew, HAR_controls_df$overall_max_skew)

HAR_skew_barplot <- ggplot(data = HAR_comp_df, aes(x = type, y = max_skew)) + 
  geom_boxplot() + labs(x = NULL, y = "|Maximum Log2 Skew|")

#check if within active HAQER/HAQER controls, there is more diff. activity in HARs
HAQER_active_df <- all_active_df %>% dplyr::filter(!is.na(HAQERs_Overlapped))%>% 
  mutate(type = "HAQER") %>% dplyr::select(region, HAQERs_Overlapped) %>% unique()

HAQER_active_controls_df <- all_active_df %>% 
  dplyr::filter(is.na(HAQERs_Overlapped)) %>%
  dplyr::filter((bp_diffs > 7) |
                  (align_type == "chimp_deletion" & deletion_size > 1)) %>% 
  mutate(type = "CONTROL") %>% dplyr::select(region, HAQERs_Overlapped) %>% unique()

HAQER_active_comp_df <- rbind(HAQER_active_df, HAQER_active_controls_df)
#add emvar information
HAQER_active_comp_df$emvar <- "Active"
HAQER_active_comp_df$emvar[which(HAQER_active_comp_df $region %in% all_emvars_df$region)] <- "Diff. Active"

HAQER_active_comp_table <- table(!is.na(HAQER_active_comp_df$HAQERs_Overlapped), HAQER_active_comp_df$emvar)
fisher.test(HAQER_active_comp_table)

HAQERs_emvars_df <- all_emvars_df %>% 
  dplyr::filter(!is.na(HAQERs_Overlapped)) %>% 
  mutate(type = "HAQER") %>% dplyr::select(region, max_skew) %>% 
  group_by(region) %>% 
  mutate(overall_max_skew = max(max_skew)) %>% 
  ungroup %>% 
  dplyr::select(region, overall_max_skew) %>%
  unique()
HAQERs_emvars_df$Overlap = "HAQER"
HAQERs_emvars_df$comp = "HAQER"

HAQERs_controls_df <- all_emvars_df %>% 
  dplyr::filter(is.na(HAQERs_Overlapped)) %>% 
  dplyr::filter((bp_diffs > 7) |
                  (align_type == "chimp_deletion" & deletion_size > 1)) %>% 
  mutate(type = "CONTROL") %>% dplyr::select(region, max_skew) %>% 
  group_by(region) %>% 
  mutate(overall_max_skew = max(max_skew)) %>% 
  ungroup %>% 
  dplyr::select(region, overall_max_skew) %>%
  unique()
HAQERs_controls_df$Overlap = "CTRL"
HAQERs_controls_df$comp = "HAQER"

HAQERs_comp_df <- rbind(HAQERs_emvars_df, HAQERs_controls_df)

wilcox.test(HAQERs_emvars_df$overall_max_skew, HAQERs_controls_df$overall_max_skew)

#make plots showing comparison to control sets
active_comp_df <- rbind(HAR_active_comp_table, HAQER_active_comp_table) %>% 
  as.data.frame()
active_comp_df$Overlap <- c("CTRL", "HAR", "CTRL", "HAQER")
active_comp_df$comp <- c("HAR", "HAR", "HAQER", "HAQER")
active_comp_df2 <- active_comp_df %>% 
  pivot_longer(col = c(Active, `Diff. Active`), 
               names_to = "activity", 
               values_to = "count")


active_comp_df2$Overlap <- factor(active_comp_df2$Overlap, levels = c("HAR", "HAQER", "CTRL"))
active_comp_df2$comp <- factor(active_comp_df2$comp, levels = c("HAR", "HAQER"))

HAR_HAQER_contig_table <- ggplot(data = active_comp_df2, 
       aes(fill = activity, y = count, x = Overlap)) + 
  geom_bar(position="fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = count), position=position_fill(vjust = 0.5), size = 2.5) +
  facet_grid(cols = vars(comp), scales = "free_x") +
  labs(x = NULL, y = "Percent of regions", fill  = "Activity") + 
  scale_fill_manual(values = c("#56B4E9", "#0079B2")) +
  theme(legend.position = "bottom")

skew_comp_df <- rbind(
  HAR_comp_df,
  HAQERs_comp_df)

#add facets
skew_comp_df$Overlap <- factor(skew_comp_df$Overlap, levels = c("HAR", "HAQER", "CTRL"))
skew_comp_df$comp <- factor(skew_comp_df$comp, levels = c("HAR", "HAQER"))

#generate plot
HAR_HAQER_skew_plot <- ggplot(data = skew_comp_df, aes(x = Overlap, y = overall_max_skew, fill = Overlap)) + 
  facet_grid(cols = vars(comp), scales = "free_x") +
  geom_boxplot(outliers = F) + 
  labs(x = NULL, y = "|Maximum Log2 Skew|", fill = NULL ) +
  annotate("text", label='NS', x = 1.5, y = 5) +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
  scale_fill_manual(values = c("darkorchid2", "goldenrod", "grey"))+
  ylim(0, 5.5)

#test if number of base pair differences between tiles predicts differential activity
bp_diffs <- data.frame(bp_diffs = c(na.omit(all_active_df$bp_diffs), na.omit(all_emvars_df$bp_diffs)), 
                       Set = c(rep(x = "Active", length(na.omit(all_active_df$bp_diffs))), rep(x = "Differentially Active", length(na.omit(all_emvars_df$bp_diffs))))
)

wilcox.test(na.omit(all_active_df$bp_diffs), na.omit(all_emvars_df$bp_diffs))
wilcox.test(na.omit(all_active_df$bp_diffs), na.omit(cartilage_metadata_experimental$bp_diffs))
wilcox.test(na.omit(all_emvars_df$bp_diffs), na.omit(cartilage_metadata_experimental$bp_diffs))

bp_diffs_activity <- ggplot(data = bp_diffs, aes(x = Set, y = bp_diffs, fill = Set))+
  geom_boxplot(outliers = F) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = "SNVs per region") + 
  scale_fill_manual(values = c("#56B4E9", "#0079B2")) + 
  annotate("text", label='***', x = 1.5, y = 9.5) + 
  ylim(0, 10)

#plot skew against bp difs
bp_diff_plot <- ggplot(data = all_emvars_df %>% 
                         dplyr::filter(align_type == "SNPs_only", skew_sign_match) %>% 
                         mutate(x_bins = cut(bp_diffs, breaks = c(0,5,10,20,45) )), 
                       aes(x = x_bins, y = max_skew)) +
  geom_boxplot(fill = "lightgrey", outliers = F) + 
  labs(x = "Number of human-chimp SNVs", y = "|Maximum log2 skew|") +
  scale_x_discrete(labels = c('1-5','5-10','11-20', '22-42')) + 
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) + annotate("text", label='***', x = 2.5, y = 5.5) + 
  ylim(0, 6)

cor.test(all_emvars_df %>% dplyr::filter(align_type == "SNPs_only", skew_sign_match) %>% pull(max_skew), 
         all_emvars_df %>% dplyr::filter(align_type == "SNPs_only", skew_sign_match) %>% pull(bp_diffs))

cor.test(all_emvars_df %>% dplyr::filter(align_type == "chimp_deletion", skew_sign_match) %>% pull(max_skew), 
         all_emvars_df %>% dplyr::filter(align_type == "chimp_deletion", skew_sign_match) %>% pull(deletion_size))

#plot skew against deletion size
deletion_size_plot <- ggplot(data = all_emvars_df %>% 
                               dplyr::filter(align_type == "chimp_deletion", skew_sign_match) %>% mutate(x_bins = cut(deletion_size, breaks = c(0,2,5,10,15,20,27))), 
                             aes(x = x_bins, y = max_skew)) + 
  geom_boxplot() +
  
  labs(x = "Size of chimp deletion(s) (base pairs)", y = "|Maximum log2 skew|") +
  scale_x_discrete(labels = c('1-2','3-5', '5-10','11-15', "15-20", "21-26")) + 
  scale_y_break(c(11, 18)) +
  theme(axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank())

#analyse variants with only a single human chimp difference

#number of tested tiles with a single SNV
length(which(cartilage_metadata_experimental$bp_diffs == 1 & cartilage_metadata_experimental$align_type == "SNPs_only",))

#number of tested tiles with a single SNV
length(which(all_emvars_df$bp_diffs == 1))

#now to look at PhyloP
find_SNP_pos <- function(x, y){
  count <- 0
  x_bps <- unlist(str_split(x, pattern = ""))
  y_bps <- unlist(str_split(y, pattern = ""))
  if(length(x_bps) == length(y_bps)){
    count <- which(sapply(seq(length(x_bps)),
                          function(i){
                            x_bps[i] != y_bps[i]
                          }))
    return(count)
  }else{
    return(999)
  }
}

get_SNP_pos <- function(x, pos){
  x_bps <- unlist(str_split(x, pattern = ""))
  return(x_bps[pos])
}

all_emvars_df %>% 
  dplyr::filter(align_type == "SNPs_only", bp_diffs == 1) %>% 
  rowwise() %>% 
  mutate(SNP_pos  = find_SNP_pos(human, chimp)) %>% 
  dplyr::select(SNP:start, human, chimp, SNP_pos) %>%
  mutate(human_SNP  = get_SNP_pos(human, SNP_pos), 
         chimp_SNP  = get_SNP_pos(chimp, SNP_pos)) %>% 
  mutate(SNP_start  = (start + SNP_pos -16)) %>% 
  mutate(SNP_end = SNP_start) %>% 
  dplyr::select(chr, SNP_start, SNP_end, SNP) %>% write_bed(filename = "~/Desktop/290k_MPRA_Cluster_Results/beds/all_1bp_emvar_SNPs_hg38.bed", ncol =4)

all_active_df %>% 
  dplyr::filter(align_type == "SNPs_only", bp_diffs == 1) %>% 
  rowwise() %>% 
  mutate(SNP_pos  = find_SNP_pos(human, chimp)) %>% 
  dplyr::select(SNP:start, human, chimp, SNP_pos) %>%
  mutate(human_SNP  = get_SNP_pos(human, SNP_pos), 
         chimp_SNP  = get_SNP_pos(chimp, SNP_pos)) %>% 
  mutate(SNP_start  = (start + SNP_pos -16)) %>% 
  mutate(SNP_end = SNP_start) %>% 
  dplyr::select(chr, SNP_start, SNP_end, SNP) %>% write_bed(filename = "~/Desktop/290k_MPRA_Cluster_Results/beds/all_1bp_active_SNPs_hg38.bed", ncol =4)

#phyloP scores were calculated on a separate cluster using the Zoonomia phyloP scores
#sort -k1,1 -k2,2n all_1bp_emvar_SNPs_hg38.bed | bedtools intersect -a stdin -b /cold/aokamoto/241-mammalian-2020v2.bed -sorted -wa -wb | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,9 -o distinct,mean > all_1bp_emvar_SNPs_phyloP_hg38.bed 
#sort -k1,1 -k2,2n all_1bp_emvar_SNPs_hg38.bed | bedtools intersect -a stdin -b /cold/aokamoto/241-mammalian-2020v2.bed -sorted -wa -wb | sort -k1,1 -k2,2n > all_1bp_emvar_SNPs_phyloP_hg38.bed 

PhyloP_SNP_raw <- read.table("~/Desktop/290k_MPRA_Cluster_Results/beds/all_1bp_emvar_SNPs_phyloP_hg38.bed") %>% 
  dplyr::filter(V7 ==V2) %>% dplyr::select(V4, V9)
colnames(PhyloP_SNP_raw) <- c("SNP", "phyloP")

PhyloP_SNP_df <-all_emvars_df %>% 
  dplyr::filter(align_type == "SNPs_only", bp_diffs == 1) %>% merge(y = PhyloP_SNP_raw)

PhyloP_SNP_plot <- ggplot(PhyloP_SNP_df, aes(x = phyloP, y = max_skew)) + 
  geom_point() + labs(x = "phyloP at human-chimp SNV", y = "|Maximum Log2 Skew|")

cor.test(PhyloP_SNP_df$phyloP, PhyloP_SNP_df$max_skew)

active_PhyloP_SNP_raw <- read.table("~/Desktop/290k_MPRA_Cluster_Results/beds/all_1bp_active_SNPs_phyloP_hg38.bed") %>% 
  dplyr::filter(V7 ==V2) %>% dplyr::select(V4, V9)
colnames(active_PhyloP_SNP_raw) <- c("SNP", "phyloP")
active_PhyloP_SNP_raw$emvar <- "Active"
active_PhyloP_SNP_raw$emvar[which(active_PhyloP_SNP_raw$SNP %in% all_emvars_df$SNP)] <- "Diff. Active"

wilcox.test(active_PhyloP_SNP_raw$phyloP[which(active_PhyloP_SNP_raw$emvar == "Active")], 
            active_PhyloP_SNP_raw$phyloP[which(active_PhyloP_SNP_raw$emvar == "Diff. Active")])

active_phyloP_plot <- ggplot(data = active_PhyloP_SNP_raw, aes(x = emvar, y = phyloP, fill = emvar))+
  geom_boxplot(outliers = F) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = "PhyloP at SNV") + 
  scale_fill_manual(values = c("#56B4E9", "#0079B2")) + 
  annotate("text", label='NS', x = 1.5, y = 2.5) + 
  ylim(-3, 3)

#motifbreakr
#to run motifbreakR, I needed to install a whole bunch of other packages
library(stringr)
library(TFMPvalue)
library(matrixStats)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)

all_emvars_df_SNVs <- all_emvars_df %>% dplyr::filter(bp_diffs ==1) %>% 
  rowwise() %>% 
  mutate(SNP_pos  = find_SNP_pos(human, chimp)) %>% 
  mutate(human_bp = unlist(str_split(human, pattern = ""))[SNP_pos], 
         chimp_bp = unlist(str_split(chimp, pattern = ""))[SNP_pos], 
         SNV_id = paste(chr, start+SNP_pos-16, human_bp, chimp_bp, sep = ":"), strand = "*", score2 = 0, start2 = start+SNP_pos-17, end2 = start+SNP_pos-16)

all_emvars_df_SNVs %>% 
  rowwise() %>% 
  dplyr::select(chr, start2, end2, SNV_id, score2, strand) %>% 
  write_bed(filename = "~/Desktop/290k_MPRA_Cluster_Results/beds/all_emvars_df_SNVs_monaLisa.bed", ncol =6)



snps.mb.frombed <- snps.from.file(file = "~/Desktop/290k_MPRA_Cluster_Results/beds/all_emvars_df_SNVs_monaLisa.bed",
                                  #dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                  search.genome = BSgenome.Hsapiens.UCSC.hg38,
                                  format = "bed", check.unnamed.for.rsid = F)

motifbreakr.results <- motifbreakR(snpList = snps.mb.frombed, pwmList = subset(subset(MotifDb, dataSource == "HOCOMOCOv10"), 
                                                                               organism == "Hsapiens"), 
                                   filterp = TRUE, 
                                   threshold = 1e-4,
                                   method = "ic",
                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                   BPPARAM = BiocParallel::SerialParam())


all_emvars_df_SNVs$strong_TFs <- NA
all_emvars_df_SNVs$weak_TFs <- NA
all_emvars_df_SNVs$all_TFs <- NA
all_emvars_df_SNVs$strong_TFs_n <- NA
all_emvars_df_SNVs$weak_TFs_n <- NA
all_emvars_df_SNVs$all_TFs_n <- NA

for (snp_pos in 1:nrow(all_emvars_df_SNVs)){
  #identify strong and weak SNPs
  all_emvars_df_SNVs$strong_TFs[snp_pos] <- paste(unique(motifbreakr.results$geneSymbol[which(motifbreakr.results$SNP_id == all_emvars_df_SNVs$SNV_id[snp_pos] & motifbreakr.results$effect == "strong")]), collapse = ",")
  all_emvars_df_SNVs$weak_TFs[snp_pos] <- paste(unique(motifbreakr.results$geneSymbol[which(motifbreakr.results$SNP_id == all_emvars_df_SNVs$SNV_id[snp_pos] & motifbreakr.results$effect == "weak")]), collapse = ",")
  all_emvars_df_SNVs$all_TFs[snp_pos] <- paste(unique(motifbreakr.results$geneSymbol[which(motifbreakr.results$SNP_id == all_emvars_df_SNVs$SNV_id[snp_pos])]), collapse = ",")
  all_emvars_df_SNVs$strong_TFs_n[snp_pos] <- length(unique(motifbreakr.results$geneSymbol[which(motifbreakr.results$SNP_id == all_emvars_df_SNVs$SNV_id[snp_pos] & motifbreakr.results$effect == "strong")]))
  all_emvars_df_SNVs$weak_TFs_n[snp_pos] <- length(unique(motifbreakr.results$geneSymbol[which(motifbreakr.results$SNP_id == all_emvars_df_SNVs$SNV_id[snp_pos] & motifbreakr.results$effect == "weak")]))
  all_emvars_df_SNVs$all_TFs_n[snp_pos]  <- length(unique(motifbreakr.results$geneSymbol[which(motifbreakr.results$SNP_id == all_emvars_df_SNVs$SNV_id[snp_pos])]))
}

strong_tfs <- unlist(str_split(all_emvars_df_SNVs$strong_TFs, ",")) %>% table() %>% as.data.frame() %>% arrange(Freq)
weak_tfs <- unlist(str_split(all_emvars_df_SNVs$weak_TFs, ",")) %>% table() %>% as.data.frame() %>% arrange(Freq)

write.table(x = strong_tfs, file = "~/Desktop/290k_MPRA_Cluster_Results/results/single_SNV_strong_TFs.txt")
write.table(x = weak_tfs, file = "~/Desktop/290k_MPRA_Cluster_Results/results/single_SNV_weak_TFs.txt")

#load and format TF data
strong_tfs <- read.table(file = "~/Desktop/290k_MPRA_Cluster_Results/results/single_SNV_strong_TFs.txt")
weak_tfs <- read.table(file = "~/Desktop/290k_MPRA_Cluster_Results/results/single_SNV_weak_TFs.txt")
colnames(strong_tfs) <- c("TF", "Count")
strong_tfs$type = "Strong"
colnames(weak_tfs) <- c("TF", "Count")
weak_tfs$type = "Weak"

tf_total <- merge(strong_tfs, weak_tfs, by = "TF", all.x =T, all.y = T)
tf_total[is.na(tf_total)] <- 0
tf_total$total <- tf_total$Count.x + tf_total$Count.y
tf_total$rank <- rank(tf_total$total,ties.method = "random")
tf_total_sub <- tf_total %>% dplyr::select(TF, total, rank)
top20_TFs <- tf_total_sub %>% dplyr::filter(TF != "") %>% arrange(desc(rank)) %>% head(n = 20) %>% pull(TF)

#merge data and plot
disrupted_TFs <- rbind(strong_tfs, weak_tfs)
disrupted_TFs_plotting <- disrupted_TFs %>% left_join(tf_total_sub) %>% 
  dplyr::filter(TF != "") %>% dplyr::filter(TF %in% top20_TFs)
disrupted_TFs_plotting$TF <- factor(disrupted_TFs_plotting$TF, levels = top20_TFs)

TF_disruption_plot <- ggplot(data = disrupted_TFs_plotting, 
       aes(x = TF, y = Count, fill = type)) + 
         geom_bar(position = "stack", stat = "identity") + 
  theme_classic() +
  scale_fill_manual(values = c("grey30", "grey")) +
  geom_text(aes(label = TF, y = total +1, angle = 90, hjust = 0)) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = c(0.8, 0.9)) + 
  labs(x = NULL, fill = "TF Disruption") +
  ylim(0,30) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

plot_h = 0.25

#create the figure
HAR_HAQER_total_plot <- ggdraw() +
  draw_plot(HAR_HAQER_contig_table , x = 0, y = plot_h*3, width = 0.5, height = plot_h) +
  draw_plot(HAR_HAQER_skew_plot, x = 0.5, y = plot_h*3, width = 0.5, height = plot_h) +
  draw_plot(bp_diffs_activity, x = 0, y = plot_h*2, width = 0.5, height = plot_h) +
  draw_plot(bp_diff_plot, x = 0.5, y = plot_h*2, width = 0.5, height = plot_h) +
  draw_plot(active_phyloP_plot, x = 0, y = plot_h, width = 0.5, height = plot_h) +
  draw_plot(PhyloP_SNP_plot, x = 0.5, y = plot_h, width = 0.5, height = plot_h) +
  draw_plot(TF_disruption_plot, x = 0, y = 0, width = 1, height = plot_h) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F", "G"), 
                  size = 10,
                  x = c(0, 0.5, 0, 0.5, 0, 0.5, 0), 
                  y = c(plot_h*4, plot_h*4, plot_h*3, plot_h*3, plot_h*2, plot_h*2, plot_h)) + 
  bgcolor("white")

ggsave(plot = HAR_HAQER_total_plot + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Cartilage MPRA Paper/Code/fig4_evol_sequence_features.png", 
       device = "png", dpi = 300, height = 9.5, width = 6.5, units = "in")