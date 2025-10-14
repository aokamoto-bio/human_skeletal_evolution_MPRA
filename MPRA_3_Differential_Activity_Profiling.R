## Code for Massively parallel functional screen identifies thousands of regulatory differences in human versus chimpanzee postcranial skeletal development

## Part 3: Characteristics of differentially active sequences

#load data
#modify to specify file path for your computer
setwd("~/Desktop/290k_MPRA_Cluster_Results/results/")
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggVennDiagram)
library(cowplot)
library(magick) # for adding images to ggplots

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

#compare emVar skews across cell lines
CHON_TC28_emvar_tiles <- CHON002_emVars$SNP[which(CHON002_emVars$SNP %in% TC28_emVars$SNP)]
CHON_K562_emvar_tiles <- CHON002_emVars$SNP[which(CHON002_emVars$SNP %in% K562_emVars$SNP)]
TC28_K562_emvar_tiles <- TC28_emVars$SNP[which(TC28_emVars$SNP %in% K562_emVars$SNP)]

#build scatterplots to show results
cellline_skew_corr <- data.frame(tile1_skew = c(CHON002_emVars$Log2Skew[which(CHON002_emVars$SNP %in% CHON_TC28_emvar_tiles)], 
                                                CHON002_emVars$Log2Skew[which(CHON002_emVars$SNP %in% CHON_K562_emvar_tiles)], 
                                                TC28_emVars$Log2Skew[which(TC28_emVars$SNP %in% TC28_K562_emvar_tiles)]), 
                                         tile2_skew = c(TC28_emVars$Log2Skew[which(TC28_emVars$SNP %in% CHON_TC28_emvar_tiles)], 
                                                        K562_emVars$Log2Skew[which(K562_emVars$SNP %in% CHON_K562_emvar_tiles)], 
                                                        K562_emVars$Log2Skew[which(K562_emVars$SNP %in% TC28_K562_emvar_tiles)]), 
                                 cellline = c(rep("CHON002 vs TC28", times = length(CHON_TC28_emvar_tiles)), 
                                              rep("CHON002 vs K562", times = length(CHON_K562_emvar_tiles)), 
                                              rep("TC28 vs K562", times = length(TC28_K562_emvar_tiles)))
)


#create the plot
cellline_skew_corr_plot <- ggplot(data = cellline_skew_corr , 
                                          mapping = aes(x = tile1_skew, y = tile2_skew)) + 
  geom_point(size = 0.2) + 
  facet_grid(~cellline) + 
  stat_cor(method="pearson", label.y = 11, p.accuracy = 0.00001) + 
  labs(x = "Cell type 1 log2 fold change skew", y = "Cell type 2 log2 fold\nchange skew") +
  ylim(-12, 12)

#count total number of regions with both tiles emVars
(n_2tile_emvar <- rbind(CHON002_emVars, K562_emVars,TC28_emVars) %>% 
  dplyr::select(region, tile) %>% 
  unique() %>% 
  group_by(region) %>% 
  summarise(n =n()) %>% 
  dplyr::filter(n > 1) %>% 
    nrow())

#count total number of regions with only one tile an emVar
(n_1tile_emvar <- rbind(CHON002_emVars, K562_emVars,TC28_emVars) %>% 
  dplyr::select(region, tile) %>% 
  unique() %>% 
  group_by(region) %>% 
  summarise(n =n()) %>% 
  dplyr::filter(n == 1) %>% 
    nrow())

#total emvars
n_2tile_emvar + n_1tile_emvar

#deeper emvar 2 tile analysis
all_2tile_emvar <- rbind(CHON002_emVars, K562_emVars, TC28_emVars) %>% 
  dplyr::select(region, tile) %>% 
  unique() %>% 
  group_by(region) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n > 1) %>% pull(region)

chon_2tile_emvar <- rbind(CHON002_emVars) %>% 
  dplyr::select(region, tile) %>% 
  unique() %>% 
  group_by(region) %>% 
  dplyr::summarise(n =n()) %>% 
  dplyr::filter(n > 1) %>% pull(region)

k562_2tile_emvar <- rbind(K562_emVars) %>% 
  dplyr::select(region, tile) %>% 
  unique() %>% 
  group_by(region) %>% 
  dplyr::summarise(n =n()) %>% 
  dplyr::filter(n > 1) %>% pull(region)

tc28_2tile_emvar <- rbind(TC28_emVars) %>% 
  dplyr::select(region, tile) %>% 
  unique() %>% 
  group_by(region) %>% 
  dplyr::summarise(n =n()) %>% 
  dplyr::filter(n > 1) %>% pull(region)


#count regions with 2 regions across two cell types
length(all_2tile_emvar)-length(unique(c(tc28_2tile_emvar, chon_2tile_emvar, k562_2tile_emvar)))
tiles_split_emvar <- all_2tile_emvar[!all_2tile_emvar %in% unique(c(tc28_2tile_emvar, chon_2tile_emvar, k562_2tile_emvar))]
length(intersect(tiles_split_emvar[which(tiles_split_emvar %in% TC28_emVars$region)], tiles_split_emvar[which(tiles_split_emvar %in% CHON002_emVars$region)]))
length(intersect(tiles_split_emvar[which(tiles_split_emvar %in% TC28_emVars$region)], tiles_split_emvar[which(tiles_split_emvar %in% K562_emVars$region)]))
length(intersect(tiles_split_emvar[which(tiles_split_emvar %in% CHON002_emVars$region)], tiles_split_emvar[which(tiles_split_emvar %in% K562_emVars$region)]))

TC28_emvar_tiles <- TC28_emVars %>% 
  dplyr::filter(region %in% tc28_2tile_emvar) %>% 
  dplyr::select(region, tile, Log2Skew) %>% distinct() %>% 
  group_by(region) %>% 
  pivot_wider(values_from = Log2Skew, names_from = tile)

CHON002_emvar_tiles <- CHON002_emVars %>% 
  dplyr::filter(region %in% chon_2tile_emvar) %>% 
  dplyr::select(region, tile, Log2Skew) %>% distinct() %>% 
  group_by(region) %>% 
  pivot_wider(values_from = Log2Skew, names_from = tile)

K562_emvar_tiles <- K562_emVars %>% 
  dplyr::filter(region %in% k562_2tile_emvar) %>% 
  dplyr::select(region, tile, Log2Skew) %>% distinct() %>% 
  group_by(region) %>% 
  pivot_wider(values_from = Log2Skew, names_from = tile)

cor.test(TC28_emvar_tiles$`1`, TC28_emvar_tiles$`2`)
cor.test(CHON002_emvar_tiles$`1`, CHON002_emvar_tiles$`2`)
cor.test(K562_emvar_tiles$`1`, K562_emvar_tiles$`2`)

#build scatterplots to show results
two_tile_simple_corr_emvar <- data.frame(tile1_skew = c(TC28_emvar_tiles$`1`, CHON002_emvar_tiles$`1`, K562_emvar_tiles$`1`), 
                                       tile2_skew = c(TC28_emvar_tiles$`2`, CHON002_emvar_tiles$`2`, K562_emvar_tiles$`2`), 
                                       cellline = c(rep("TC28", times = nrow(TC28_emvar_tiles)), 
                                                    rep("CHON002", times = nrow(CHON002_emvar_tiles)), 
                                                    rep("K562", times = nrow(K562_emvar_tiles)))
)


#create the plot
two_tile_simple_emvar_corr_plot <- ggplot(data = two_tile_simple_corr_emvar, 
                                    mapping = aes(x = tile1_skew, y = tile2_skew)) + 
  geom_point(size = 0.2) + 
  facet_grid(~cellline) + 
  stat_cor(method="pearson", label.y = 10) + 
  labs(x = "Tile 1 log2 fold change skew", y = "Tile 2 log2 fold\nchange skew") +
  ylim(c(-7, 11))

TC28_emvar_2tiles_sign_match <- TC28_emvar_tiles$region[which(sign(TC28_emvar_tiles$`1`) ==sign(TC28_emvar_tiles$`2`))]
CHON002_emvar_2tiles_sign_match <- CHON002_emvar_tiles$region[which(sign(CHON002_emvar_tiles$`1`) ==sign(CHON002_emvar_tiles$`2`))]
K562_emvar_2tiles_sign_match <- K562_emvar_tiles$region[which(sign(K562_emvar_tiles$`1`) ==sign(K562_emvar_tiles$`2`))]

#now for the fancier case
#compare pairwise sets 
split_tiles_cell_type_emvar <- function(tiles_split_list, cell_df1, cell_df2){
  #get list of split tiles in these two dfs
  split_tiles_emvar <- intersect(tiles_split_list[which(tiles_split_list %in% cell_df1$region)], tiles_split_list[which(tiles_split_list %in% cell_df2$region)])
  
  split_tiles_df_emvar <- data.frame(SNP = cartilage_metadata_experimental$SNP[which(cartilage_metadata_experimental$region %in% split_tiles_emvar)],
                                     region = cartilage_metadata_experimental$region[which(cartilage_metadata_experimental$region %in% split_tiles_emvar)], 
                                     tile = cartilage_metadata_experimental$tile[which(cartilage_metadata_experimental$region %in% split_tiles_emvar)], 
                                     Log2Skew = NA)
  
  #add emvar tile info from relevant cell type
  split_tiles_df_emvar$Log2Skew[which(split_tiles_df_emvar$SNP %in% cell_df1$SNP)] <- cell_df1$Log2Skew[which(cell_df1$SNP %in% split_tiles_df_emvar$SNP)]
  
  split_tiles_df_emvar$Log2Skew[which(split_tiles_df_emvar$SNP %in% cell_df2$SNP)] <- cell_df2$Log2Skew[which(cell_df2$SNP %in% split_tiles_df_emvar$SNP)]
  
  split_tiles_df_emvar2 <- split_tiles_df_emvar %>% 
    dplyr::select(-SNP) %>% 
    dplyr::filter(region %in% split_tiles_emvar) %>% 
    dplyr::select(region, tile, Log2Skew) %>% 
    distinct() %>% 
    group_by(region) %>% 
    pivot_wider(values_from = Log2Skew, names_from = tile) %>% 
    drop_na()
  
  #return filled in dataframe
  return(split_tiles_df_emvar2)
}

#chon versus tc28
split_tiles_df_chon_tc28_emvar <- split_tiles_cell_type_emvar(tiles_split_list = tiles_split_emvar, cell_df1 = TC28_emVars, cell_df2 = CHON002_emVars)
split_tiles_df_chon_k562_emvar <- split_tiles_cell_type_emvar(tiles_split_list = tiles_split_emvar, cell_df1 = K562_emVars, cell_df2 = CHON002_emVars)
split_tiles_df_tc28_k562_emvar <- split_tiles_cell_type_emvar(tiles_split_list = tiles_split_emvar, cell_df1 = K562_emVars, cell_df2 = TC28_emVars)
length(unique(c(split_tiles_df_chon_tc28_emvar$region, split_tiles_df_chon_k562_emvar$region, split_tiles_df_tc28_k562_emvar$region)))

#build scatterplots to show results
  two_tile_complex_corr_emvar <- data.frame(tile1_skew = c(split_tiles_df_chon_tc28_emvar$`1`, split_tiles_df_chon_k562_emvar$`1`, split_tiles_df_tc28_k562_emvar$`1`), 
                                        tile2_skew = c(split_tiles_df_chon_tc28_emvar$`2`, split_tiles_df_chon_k562_emvar$`2`, split_tiles_df_tc28_k562_emvar$`2`), 
                                        cellline = c(rep("CHON002 vs TC28", times = nrow(split_tiles_df_chon_tc28_emvar)), 
                                                     rep("CHON002 vs K562", times = nrow(split_tiles_df_chon_k562_emvar)), 
                                                     rep("TC28 vs K562", times = nrow(split_tiles_df_tc28_k562_emvar)))
)

#create the plot
two_tile_complex_emvar_corr_plot <- ggplot(data = two_tile_complex_corr_emvar, 
                                     mapping = aes(x = tile1_skew, y = tile2_skew)) + 
  geom_point(size = 0.2) + 
  facet_wrap(~cellline, scales= "free") + 
  stat_cor(method="pearson", label.y = 13) + 
  labs(x = "Cell type 1 log2 fold change skew", y = "Cell type 2 log2 fold\nchange skew") +
  ylim(c(-7, 14))

#show skew
emvar_skew_bias <-  data.frame(skew = c(TC28_emVars$Log2Skew, CHON002_emVars$Log2Skew, K562_emVars$Log2Skew), 
                               cellline = c(rep("TC28", times = nrow(TC28_emVars)), 
                                            rep("CHON002", times = nrow(CHON002_emVars)), 
                                            rep("K562", times = nrow(K562_emVars)))
)
emvar_skew_bias$`Species bias` <- "Human"
emvar_skew_bias$`Species bias`[which(emvar_skew_bias$skew > 0 )] <- "Chimp"

emvar_skew_bias_plot <- ggplot(emvar_skew_bias, aes(x = skew, fill =  `Species bias`)) + 
  geom_histogram(binwidth = 0.5) + 
  facet_wrap(~cellline, scales = "free") + 
  labs(y = "Count", x = "log2 fold change skew") + 
  scale_fill_manual(, values = c("#117733","#0079B2"))

#create the supplementary figure

two_tile_figure_emvar <- ggdraw() +
  #A, tiles diff active between cell types
  draw_plot(cellline_skew_corr_plot, x = 0, y = 0.75, width = 1, height = 0.25) +
  #B, histogram tile skew
  draw_plot(emvar_skew_bias_plot, x = 0, y = 0.5, width = 1, height = 0.25) +
  #C, differential activity in both tiles within region within cell type
  draw_plot(two_tile_simple_emvar_corr_plot, x = 0, y = 0.25, width = 1, height = 0.25) +
  #C, differential activity in both tiles within region across cell types
  draw_plot(two_tile_complex_emvar_corr_plot, x = 0, y = 0, width = 1, height = 0.25) +
  #add labels
  draw_plot_label(label = c("A", "B", "C", "D"), 
                  size = 12,
                  x = c(0, 0, 0, 0), 
                  y = c(1, 0.75, 0.5, 0.25)) + 
  bgcolor("white")

ggsave(plot = two_tile_figure_emvar, 
       filename = "~/Dropbox/Cartilage MPRA Paper/Code/S7_two_tile_dif_active_scatterplots.png", 
       device = "png", dpi = 300, height = 8, width = 6, units = "in")

#plot pie chart of active sequences

#format data
emvar_pie_data <- data.frame(
  class = c("Active", "Diff. Active (2 tile)", "Diff. Active (1 tile)"), 
  count = c(active_regions_count - n_2tile_emvar - n_1tile_emvar, n_2tile_emvar, n_1tile_emvar)
)

#format for plotting
emvar_pie_data2 <- emvar_pie_data %>% 
  mutate(pct = count/sum(count)*100)
emvar_pie_data2$class <- factor(emvar_pie_data2$class, c("Diff. Active (1 tile)", "Active", "Diff. Active (2 tile)"))
emvar_pie_data2$label = c("19,212\n(62.5%)","615\n(2.0%)",  "10,909\n(35.5%)")

#generate plot
emvar_pie <- ggplot(data = emvar_pie_data2, aes(x="", y = pct, fill = class)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  #geom_label_repel(aes(label = count, fill = NULL),
  #                 size = 1.5, nudge_x = 1, show.legend = FALSE)+
  geom_text(aes(label = label, x = 1.1), position = position_stack(vjust=0.5), size = 3) +
  theme_void() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#0079B2",  "#88CCEE", "dodgerblue4"), name = NULL) + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

#check for correlation
cor.test(split_tiles_df_chon_tc28_emvar$`1`, split_tiles_df_chon_tc28_emvar$`2`)
cor.test(split_tiles_df_chon_k562_emvar$`1`, split_tiles_df_chon_k562_emvar$`2`)
cor.test(split_tiles_df_tc28_k562_emvar$`1`, split_tiles_df_tc28_k562_emvar$`2`)

ggplot(split_tiles_df_tc28_k562,
       aes(x = A_log2FC_1, y = A_log2FC_2)) +
  geom_point()

#for regions
length(unique(c(CHON002_emVars$region, K562_emVars$region, TC28_emVars$region)))

#for tiles
length(unique(c(CHON002_emVars$SNP, K562_emVars$SNP, TC28_emVars$SNP)))

# make Venn Diagrams
emVar_venn <- ggVennDiagram(
  x = list(unique(CHON002_emVars$region), unique(K562_emVars$region), unique(TC28_emVars$region)),
  category.names = c("        CHON002" , "K562" , "TC28"), 
  label_alpha = 0, label_size = 3
) + 
  labs(fill = "Count") + 
  theme(plot.title=element_text(hjust=0.5)) + 
  scale_fill_gradient(low="grey90",high = "red")

#test for bias towards human or chimp differential activity
binom.test(length(which(CHON002_emVars$Log2Skew > 0 )), n = length(CHON002_emVars$Log2Skew), p = 0.5)
binom.test(length(which(TC28_emVars$Log2Skew > 0 )), n = length(TC28_emVars$Log2Skew), p = 0.5)
binom.test(length(which(K562_emVars$Log2Skew > 0 )), n = length(K562_emVars$Log2Skew), p = 0.5)

#this requires Part 2 results be in memory

total_emvar_plot <- ggplot(data = overlap_sharing %>% dplyr::filter(celltype == "ALL"), 
                           aes(x =total_overlaps, y = emvar_overlaps)) +
  geom_point() + 
  theme_bw() + 
  labs(x = "Total regions \nper tissue", 
       y = "Regions with diff. \nactivity per tissue") + 
  stat_cor(method="pearson", p.accuracy = 1e-7, label.y = 17000) +
  theme(plot.margin=unit(c(0.5,0.5,0.5, 0.5), "cm")) +
  ylim(0, 17500)

unique_emvar_plot <- ggplot(data = overlap_sharing %>% dplyr::filter(celltype == "ALL"), 
                            aes(x = unique_overlaps, y = emvar_unique_overlaps))+
  geom_point() + 
  theme_bw() + 
  labs(x = "Unique regions \nper tissue", 
       y = "Unique regions with \ndiff. activity per tissue") + 
  stat_cor(method="pearson", p.accuracy = 1e-7, label.y = 700) +
  theme(plot.margin=unit(c(0.5,0.5,0.5, 0.5), "cm"))+
  ylim(0, 720)

emvar_regions = unique(c(CHON002_emVars$region, TC28_emVars$region, K562_emVars$region)) 
emvar_snps = unique(c(CHON002_emVars$SNP, TC28_emVars$SNP, K562_emVars$SNP)) 
cartilage_metadata_experimental$diff_activity <- "no_dif"
cartilage_metadata_experimental$diff_activity[which(cartilage_metadata_experimental$SNP %in% emvar_snps)] <- "dif"

# har+ table with percent active, emvars, avg accessibility, avg. effect size

CHON002_glm <- read.table("290k_CHON002_emVAR_glm_20250311.out", header = T)
TC28_glm <- read.table("290k_TC28_emVAR_glm_20250311.out", header = T) 
K562_glm <- read.table("290k_K562_emVAR_glm_20250311.out", header = T)

emvar_snps = unique(c(CHON002_emVars$SNP, TC28_emVars$SNP, K562_emVars$SNP)) 
active_snps = unique(c(CHON002_active$SNP, TC28_active$SNP, K562_active$SNP))

emvar_df <- cartilage_metadata_experimental %>% dplyr::filter(SNP %in% emvar_snps)
emvar_df$status <- "none"
emvar_df$status[which(!is.na(emvar_df$HARs_Overlapped))] <- "HAR"
emvar_df$status[which(!is.na(emvar_df$HAQERs_Overlapped))] <- "HAQER"
emvar_df$status[which(!is.na(emvar_df$hCONDELs_Overlapped))] <- "hCONDEL"
emvar_df$status[which(!is.na(emvar_df$hCONDELs_Overlapped)&!is.na(emvar_df$HARs_Overlapped))] <- "HAR/hCONDEL"
emvar_df$status[which(!is.na(emvar_df$HAQERs_Overlapped)&!is.na(emvar_df$HARs_Overlapped))] <- "HAR/HAQER"
emvar_df$status[which(!is.na(emvar_df$HAQERs_Overlapped)&!is.na(emvar_df$hCONDELs_Overlapped))] <- "hCONDEL/HAQER"
emvar_df$status[which(!is.na(emvar_df$HAQERs_Overlapped)&!is.na(emvar_df$HARs_Overlapped)&!is.na(emvar_df$hCONDELs_Overlapped))] <- "All_Three"

#get the most significant skew from across the three cell types
emvar_df$Skew_logFDR_act <- NA
emvar_df$Log2Skew <- NA
emvar_df$most_sig_cell <- NA

for (s in emvar_df$SNP){
  chon_log_p <- CHON002_glm$Skew_logFDR_act[which(CHON002_glm$SNP == s)] %>% replace(is.na(.), 0)
  tc28_log_p <- TC28_glm$Skew_logFDR_act[which(TC28_glm$SNP == s)] %>% replace(is.na(.), 0)
  k562_log_p <- K562_glm$Skew_logFDR_act[which(K562_glm$SNP == s)] %>% replace(is.na(.), 0)
  
  #if CHON most significant
  if((chon_log_p  > tc28_log_p) && (chon_log_p  > k562_log_p)){
    emvar_df$Skew_logFDR_act[which(emvar_df$SNP == s )] <- CHON002_glm$Skew_logFDR_act[which(CHON002_glm$SNP == s )]
    emvar_df$Log2Skew[which(emvar_df$SNP == s )] <- CHON002_glm$Log2Skew[which(CHON002_glm$SNP == s )]
    emvar_df$most_sig_cell[which(emvar_df$SNP == s )] <- "CHON002"
    
  }
  #if TC28 most significant
  if((tc28_log_p  > chon_log_p) && (tc28_log_p > k562_log_p)){
    emvar_df$Skew_logFDR_act[which(emvar_df$SNP == s )] <- TC28_glm$Skew_logFDR_act[which(TC28_glm$SNP == s )]
    emvar_df$Log2Skew[which(emvar_df$SNP == s )] <- TC28_glm$Log2Skew[which(TC28_glm$SNP == s )]
    emvar_df$most_sig_cell[which(emvar_df$SNP == s )] <- "TC28"
    
  }
  #if K562 most significant
  if((k562_log_p  > chon_log_p) && (k562_log_p > tc28_log_p)){
    emvar_df$Skew_logFDR_act[which(emvar_df$SNP == s )] <- K562_glm$Skew_logFDR_act[which(K562_glm$SNP == s )]
    emvar_df$Log2Skew[which(emvar_df$SNP == s )] <- K562_glm$Log2Skew[which(K562_glm$SNP == s )]
    emvar_df$most_sig_cell[which(emvar_df$SNP == s )] <- "K562"
    
  }
  
  
}

#%>% dplyr::filter(chimp_multi_match == "NO")
#%>% dplyr::filter(bp_diffs < 100)
ggplot(data = emvar_df %>% dplyr::filter(bp_diffs < 100), 
       aes(y = Log2Skew, x = bp_diffs, color = status)) + 
  geom_jitter(size = 0.5) 


fisher.test(table(K562_active$Log2Skew > 0, !is.na(K562_active$HAQERs_Overlapped)))


#create a function to fill out the dataframe per cell type

#return 9 values, emvars, active only, inactive, for each timepoint
get_pie_tp_data <- function(active_df, emvars_df){
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
  
  #get active 
  early_emVars <- unique(emvars_df$region[which(emvars_df$timepoint =="EARLY_ONLY")])
  late_emVars <- unique(emvars_df$region[which(emvars_df$timepoint =="LATE_ONLY")])
  shared_emVars <- unique(emvars_df$region[which(emvars_df$timepoint == "SHARED")])
  
  return(c(length(early_emVars), 
           length(early_active) - length(early_emVars),
           early_total -length(early_active),
           length(late_emVars), 
           length(late_active) - length(late_emVars),
           late_total -length(late_active),
           length(shared_emVars), 
           length(shared_active) - length(shared_emVars),
           shared_total - length(shared_active)))
}

human_tp_df$count[1:9] <- get_pie_tp_data(rbind(CHON002_active, K562_active, TC28_active), 
                                          emvars_df = rbind(CHON002_emVars, K562_emVars, TC28_emVars))
human_tp_df$count[10:18] <- get_pie_tp_data(CHON002_active, emvars_df = CHON002_emVars)
human_tp_df$count[19:27] <- get_pie_tp_data(K562_active, emvars_df = K562_emVars)
human_tp_df$count[28:36]<- get_pie_tp_data(TC28_active, emvars_df = TC28_emVars)

human_tp_df$class[which(human_tp_df$class == "emVar")] <- "Diff. Active"


human_tp_celltype_pie <- ggplot(data = human_tp_df %>% 
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
  scale_fill_manual(values = c("#0079B2", "#56B4E9", "grey"), name = NULL)

ggsave(plot = human_tp_celltype_pie, filename = "~/Desktop/290k_MPRA_Cluster_Results/results/290k_MPRA_human_timepoint_celltype_pie.png", device = "png", dpi = 300, height = 5, width = 5, units = "in", bg="white")



#compare emvars
#now direct species comparison
#repeat for chimp
CHON_TFs_human_chimp_emvar <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% CHON002_emVars$SNP)], 
                                                        sequences2 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% CHON002_emVars$SNP)], 
                                                        ids = c("CHON_diff_active_human", "CHON_diff_active_chimp")
) %>% as_tibble() %>% dplyr::filter(log2enr > 0 & negLog10Padj > 4)


TC28_TFs_human_chimp_emvar <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% TC28_emVars$SNP)], 
                                                        sequences2 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% TC28_emVars$SNP)], 
                                                        ids = c("TC28_diff_active_human", "TC28_diff_active_chimp")
) %>% as_tibble() %>% dplyr::filter(log2enr > 0 & negLog10Padj > 4)

K562_TFs_human_chimp_emvar <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% K562_emVars$SNP)], 
                                                        sequences2 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% K562_emVars$SNP)], 
                                                        ids = c("K562_diff_active_human", "K562_diff_active_chimp")
) %>% as_tibble() %>% 
  dplyr::filter(log2enr > 0 & negLog10Padj > 4)

#get frequency of number of tissues with accessibility per region
all_active_tissue_freq <- cartilage_metadata_experimental %>% 
  dplyr::filter(region %in% unique(c(K562_active$region, TC28_active$region, CHON002_active$region))) %>% 
  dplyr::select(tissue_overlap_n) %>% 
  group_by(tissue_overlap_n) %>% 
  summarize(n = n()) %>% 
  mutate(celltype = "Active")

#get frequency of number of tissues with accessibility per region
all_emvar_tissue_freq <- cartilage_metadata_experimental %>% 
  dplyr::filter(region %in% unique(c(K562_emVars$region, TC28_emVars$region, CHON002_emVars$region))) %>% 
  dplyr::select(tissue_overlap_n) %>% 
  group_by(tissue_overlap_n) %>% 
  summarize(n = n()) %>% 
  mutate(celltype = "Diff. Active")

#generate the plot
total_tissue_emvar_freq_plot <- ggplot(data = all_active_tissue_freq,
                                          aes(x = tissue_overlap_n, y = n, fill = celltype)) + 
  geom_bar(stat="identity") + 
  geom_bar(data = all_emvar_tissue_freq, stat="identity") + 
  theme_bw() + 
  labs(x = "Number of tissues with accessibility", y  = "Frequency", fill = "") + 
  theme(legend.position = c(0.5, 0.8), 
        legend.direction="horizontal", 
        legend.background = element_rect(fill="lightgray", linewidth = 0.5, linetype="solid")) + 
  scale_fill_manual(breaks=c('Diff. Active', 'Active'), values = c( "#56B4E9"), name = NULL)

human_tp_df$class <- factor(human_tp_df$class, levels = c("Differentially Active", "Active", "Inactive"))
human_tp_total_pie_emvar <- ggplot(data = human_tp_df %>% 
                                     dplyr::filter(class != "Inactive") %>% 
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
  scale_fill_manual(values = c("#0079B2", "#56B4E9"), labels = c("Diff. Active", "Active"), name = NULL) 

#phyloP for active
emvar_phylop <- read.table(file = "~/Desktop/290k_MPRA_Cluster_Results/beds/290k_active_tiles_phyloP_hg38.bed", header = F)
emvar_phylop$Set <- "Active"
colnames(emvar_phylop)[1:5] <- c("chr", "start", "stop", "SNP", "PhyloP")
emvar_phylop$Set[which(emvar_phylop$SNP %in% unique(c(CHON002_emVars$SNP, K562_emVars$SNP, TC28_emVars$SNP)))] <- "Diff. Active"
emvar_phylop$Set <- factor(emvar_phylop$Set, levels = c('Diff. Active', 'Active'))

wilcox.test(emvar_phylop$PhyloP[which(emvar_phylop$Set == "Active")], emvar_phylop$PhyloP[which(emvar_phylop$Set == "Diff. Active")])

emvar_phyloP <- ggplot(data = emvar_phylop, aes(x = Set, y = PhyloP, fill = Set))+
  #geom_violin() +
  geom_boxplot(outliers = F)  + labs(x = NULL, y = "PhyloP") + 
  scale_fill_manual(breaks=c('Diff. Active', 'Active'), values = c("#0079B2", "#56B4E9")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  annotate("text", label='***', x = 1.5, y = 1.9) + 
  ylim(-1.5, 2)


#TSS for emvars
emvar_TSS <- read.table(file = "~/Desktop/290k_MPRA_Cluster_Results/beds/290k_active_regions_TSS_hg38.bed", header = F)
emvar_TSS$Set <- "Active"

colnames(emvar_TSS)[c(1:3, 10)] <- c("chr", "start", "stop", "TSS")
emvar_TSS$Region <- paste(emvar_TSS$chr, emvar_TSS$start, sep = ":")

emvar_TSS$Set[which(emvar_TSS$Region %in% unique(c(CHON002_emVars$region, K562_emVars$region, TC28_emVars$region)))] <- "Diff. Active"
emvar_TSS$Set <- factor(emvar_TSS$Set, levels = c('Diff. Active', 'Active'))

wilcox.test(emvar_TSS$TSS[which(emvar_TSS$Set == "Active")], emvar_TSS$TSS[which(emvar_TSS$Set == "Diff. Active")])
emvar_TSS$kTSS <- emvar_TSS$TSS/1000

emvar_TSS_plot <- ggplot(data = emvar_TSS, aes(x = Set, y = kTSS, fill = Set))+
  #geom_violin() +
  geom_boxplot(outliers = F)  + 
  labs(x = NULL, y = "Distance to TSS\n(Kilobases)") + 
  scale_fill_manual(breaks=c('Diff. Active', 'Active'), values = c("#0079B2", "#56B4E9")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  annotate("text", label='*', x = 1.5, y = 30) + 
  ylim(0, 30)


#create the figure
emvar_total_plot <- ggdraw() +
  draw_plot(emvar_pie, x = 0, y = .7, width = 0.4, height = .3) +
  draw_plot(emVar_venn, x = 0.4, y = .7, width = 0.6, height = .3) +
  draw_plot(total_emvar_plot, x = 0, y = .45, width = 0.5, height = .25) +
  draw_plot(unique_emvar_plot, x = 0.5, y = .45, width = 0.5, height = .25) +
  draw_plot(total_tissue_emvar_freq_plot, x = 0, y = 0.2, width = 0.5, height = .25) +
  draw_plot(emvar_phyloP, x = 0.5, y = 0.2, width = 0.25, height = .25) +
  draw_plot(emvar_TSS_plot, x = 0.75, y = 0.2, width = 0.25, height = .25) +
  draw_plot(human_tp_total_pie_emvar, x = 0, y = 0, width = 1, height = .2) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F", "G", "H"), 
                  size = 10,
                  x = c(0, 0.4, 0, 0.5, 0, 0.5, 0.75, 0), 
                  y = c(1, 1, 0.7, 0.7, 0.45, 0.45, 0.45, 0.2)) + 
  bgcolor("white")


ggsave(plot = emvar_total_plot + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Cartilage MPRA Paper/Code/fig3_diff_activity_plot.png", 
       device = "png", dpi = 300, height = 9.5, width = 6.5, units = "in")


#finally, check TF binding differences

CHON_emvar_TFs <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% CHON002_emVars$SNP)], 
                                      sequences2 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% CHON002_active$SNP)], 
                                      ids = c("CHON_emvar", "CHON_all") )%>% 
  as_tibble() %>% dplyr::filter(negLog10Padj > 4)


TC28_emvar_TFs <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% TC28_emVars$SNP)], 
                                      sequences2 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% TC28_active$SNP)], 
                                      ids = c("TC28_emvar", "TC28_active")
) %>% as_tibble() %>% dplyr::filter(negLog10Padj > 4)

K562_emvar_TFs <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% K562_emVars$SNP)], 
                                      sequences2 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% K562_active$SNP)], 
                                      ids = c("K562_emvar", "K562_active")
) %>% as_tibble() %>% dplyr::filter(negLog10Padj > 4)

#depleted_TFs <- intersect(CHON_emvar_TFs$motif.name[which(CHON_TFs$SampleID == "all")], intersect(TC28_TFs$motif.name[which(TC28_TFs$SampleID == "all")], K562_TFs$motif.name[which(K562_TFs$SampleID == "all")]))

# make Venn Diagrams
emvar_human_TF_venn <- ggVennDiagram(
  x = list(unique(CHON_emvar_TFs$motif.name), unique(K562_emvar_TFs$motif.name), unique(TC28_emvar_TFs$motif.name)),
  category.names = c("  CHON002" , "K562" , "TC28") 
) + labs(fill = "TF Count")

#get unique TFs
chon_emvar_unique_TFs <- CHON_emvar_TFs %>%  dplyr::filter(motif.name %in% setdiff(CHON_emvar_TFs$motif.name, c(TC28_emvar_TFs$motif.name, K562_emvar_TFs$motif.name)))
tc28_emvar_unique_TFs <- TC28_emvar_TFs %>%  dplyr::filter(motif.name %in% setdiff(TC28_emvar_TFs$motif.name, c(CHON_emvar_TFs$motif.name, K562_emvar_TFs$motif.name)))
k562_emvar_unique_TFs <- K562_emvar_TFs %>%  dplyr::filter(motif.name %in% setdiff(K562_emvar_TFs$motif.name, c(TC28_emvar_TFs$motif.name, CHON_emvar_TFs$motif.name)))


#repeat for chimp
CHON_emvar_TFs_chimp <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% CHON002_emVars$SNP)], 
                                            sequences2 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% CHON002_active$SNP)], 
                                            ids = c("CHON_emvar", "CHON_active")
) %>% 
  as_tibble() %>% 
  dplyr::filter(log2enr > 0 & negLog10Padj > 4 & SampleID == "CHON_emvar")


TC28_emvar_TFs_chimp <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% TC28_emVars$SNP)], 
                                            sequences2 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% TC28_active$SNP)], 
                                            ids = c("TC28_emvar", "TC28_active")
) %>% 
  as_tibble() %>% 
  dplyr::filter(log2enr > 0 & negLog10Padj > 4 & SampleID == "TC28_emvar")

K562_emvar_TFs_chimp <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% K562_emVars$SNP)], 
                                            sequences2 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% K562_active$SNP)], 
                                            ids = c("K562_active", "all")
) %>% as_tibble() %>% dplyr::filter(log2enr > 0 & negLog10Padj > 4 & SampleID == "K562_emvar")

# make Venn Diagrams
emvar_chimp_TF_venn <- ggVennDiagram(
  x = list(unique(CHON_emvar_TFs_chimp$motif.name), unique(K562_emvar_TFs_chimp$motif.name), unique(TC28_emvar_TFs_chimp$motif.name)),
  category.names = c("  CHON002" , "K562" , "TC28") 
) + labs(fill = "TF Count")

#get unique TFs
chon_active_unique_TFs_chimp <- CHON_emvar_TFs_chimp %>%  dplyr::filter(motif.name %in% setdiff(CHON_emvar_TFs_chimp$motif.name, c(TC28_TFs_chimp$motif.name, K562_TFs_chimp$motif.name)))
tc28_active_unique_TFs_chimp <- TC28_emvar_TFs_chimp %>%  dplyr::filter(motif.name %in% setdiff(TC28_TFs_chimp$motif.name, c(CHON_TFs_chimp$motif.name, K562_TFs_chimp$motif.name)))
k562_active_unique_TFs_chimp <- K562_emvar_TFs_chimp %>%  dplyr::filter(motif.name %in% setdiff(K562_TFs_chimp$motif.name, c(TC28_TFs_chimp$motif.name, CHON_TFs_chimp$motif.name)))

#create active TF tables for supplement
TC28_emvar_TFs$celltype <- NULL
CHON_emvar_TFs$celltype <- NULL
K562_emvar_TFs$celltype <- NULL

human_emvar_TFs <- rbind(CHON_emvar_TFs, TC28_emvar_TFs, K562_emvar_TFs) %>% dplyr::select(-motif.pfm, -motif.pwm)
write_delim(human_emvar_TFs, delim = "\t", file = "/Users/alexanderokamoto/Desktop/290k_MPRA_Cluster_Results/results/human_diff_active_TFs.txt")

chimp_emvar_TFs <- rbind(CHON_emvar_TFs_chimp, TC28_emvar_TFs_chimp, K562_emvar_TFs_chimp) %>% dplyr::select(-motif.pfm, -motif.pwm)
write_delim(chimp_emvar_TFs, delim = "\t", file = "/Users/alexanderokamoto/Desktop/290k_MPRA_Cluster_Results/results/chimp_diff_active_TFs.txt")

#now direct species comparison
#repeat for chimp
CHON_emvar_TFs_human_chimp <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% CHON002_emVars$SNP)], 
                                                  sequences2 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% CHON002_emVars$SNP)], 
                                                  ids = c("CHON_emvar_human", "CHON_emvar_chimp")
) %>% as_tibble() %>% dplyr::filter(log2enr > 0 & negLog10Padj > 4)


TC28_emvar_TFs_human_chimp <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% TC28_emVars$SNP)], 
                                                  sequences2 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% TC28_emVars$SNP)], 
                                                  ids = c("TC28_emvar_human", "TC28_emvar_chimp")
) %>% as_tibble() %>% dplyr::filter(log2enr > 0 & negLog10Padj > 4)

K562_emvar_TFs_human_chimp <- pairwise_TF_analysis_MPRA(sequences1 = cartilage_metadata_experimental$human[which(cartilage_metadata_experimental$SNP %in% K562_emVars$SNP)], 
                                                  sequences2 = cartilage_metadata_experimental$chimp[which(cartilage_metadata_experimental$SNP %in% K562_emVars$SNP)], 
                                                  ids = c("K562_emvar_human", "K562_emvar_chimp")
) %>% as_tibble() %>% dplyr::filter(log2enr > 0 & negLog10Padj > 4)

# make Venn Diagrams
emvar_chimp_TF_venn <- ggVennDiagram(
  x = list(unique(CHON_TFs_chimp$motif.name), unique(K562_TFs_chimp$motif.name), unique(TC28_TFs_chimp$motif.name)),
  category.names = c("  CHON002" , "K562" , "TC28") 
) + labs(fill = "TF Count")