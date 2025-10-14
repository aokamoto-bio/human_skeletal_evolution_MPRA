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

#chimp deletions
cartilage_metadata_experimental %>% as.data.frame() %>% 
  dplyr::filter(align_type == "chimp_deletion") %>% 
  mutate(region_end = (region_start+270)) %>% 
  dplyr::select(chr, start, region_end, hCONDELs_Overlapped) %>% 
  write_bed(file = "~/Desktop/290k_MPRA_Cluster_Results/beds/290k_chimp_deletion_tiles_hg38.bed", ncol = 4)

#19,983 chimp gap regions (padded relative to human)

hcondel_phylop <- read.table("~/Desktop/290k_MPRA_Cluster_Results/beds/290k_hCONDEL_tiles_phyloP_hg38.bed") 
colnames(hcondel_phylop) <- c("chr", "start", "end", "hCONDEL", "phyloP")
hcondel_phylop$region <- paste(hcondel_phylop$chr, hcondel_phylop$start, sep = ":")

chimp_loss_phylop <- read.table("~/Desktop/290k_MPRA_Cluster_Results/beds/290k_chimp_deletion_tiles_phyloP_hg38.bed") 
colnames(chimp_loss_phylop) <- c("chr", "start2", "end", "hCONDEL", "phyloP")
#chimp_loss_phylop$region <- paste(chimp_loss_phylop$chr, chimp_loss_phylop$start2, sep = ":")
chimp_loss_phylop$length <- chimp_loss_phylop$end - chimp_loss_phylop$start2
chimp_loss_phylop2 <- chimp_loss_phylop %>% dplyr::filter(length == 270)

length(which(chimp_loss_phylop$phyloP > 0))

cartilage_metadata_experimental$start2 <- cartilage_metadata_experimental$start
cartilage_metadata_experimental$start2[which(cartilage_metadata_experimental$tile == 2 )] <- cartilage_metadata_experimental$start[which(cartilage_metadata_experimental$tile == 2 )] + 270
chimp_loss_df <- merge(cartilage_metadata_experimental, chimp_loss_phylop2, by  = c("chr", "start2"), all.y = T) %>% unique()#one annoying duplicated row

chimp_loss_df_reduced <- chimp_loss_df %>% 
  dplyr::select(SNP, deletion_size, phyloP, activity, max_skew) %>%  
  dplyr::filter(activity != "Inactive") %>% 
  dplyr::filter(phyloP > 0) %>% unique()

chimp_loss_df_reduced$activity[which(chimp_loss_df_reduced$SNP %in% all_emvars_df$SNP)] <- "Differentially Active"

wilcox.test(chimp_loss_df_reduced$deletion_size[which(chimp_loss_df_reduced$activity == "Active")], 
            chimp_loss_df_reduced$deletion_size[which(chimp_loss_df_reduced$activity == "Differentially Active")])

chimp_loss_activity_conserved <- ggplot(data = chimp_loss_df_reduced, aes(x = activity, y = deletion_size, fill = activity))+
  geom_boxplot(outliers = F) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = "Size of chimp deletion(s) (bp)") + 
  scale_fill_manual(values = c("#56B4E9", "#0079B2")) + 
  annotate("text", label='NS', x = 1.5, y = 10) + 
  ylim(0, 11)

cor.test(all_emvars_df %>% dplyr::filter(SNP %in% chimp_loss_df$SNP) %>% pull(max_skew), 
         all_emvars_df %>% dplyr::filter(SNP %in% chimp_loss_df$SNP) %>% pull(deletion_size))

#plot skew against deletion size
deletion_size_plot_conserved <- ggplot(data = all_emvars_df %>% 
                               dplyr::filter(SNP %in% chimp_loss_df$SNP) %>% 
                               mutate(x_bins = cut(deletion_size, breaks = c(0,2,5,10,15,20,27))), 
                             aes(x = x_bins, y = max_skew)) + 
  geom_boxplot(fill = "lightgrey", outliers = F) +
  labs(x = "Size of chimp deletion(s) (base pairs)", y = "|Maximum log2 skew|") +
  scale_x_discrete(labels = c('1-2','3-5', '5-10','11-15', "15-20", "21-26")) + 
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) + 
  annotate("text", label='NS', x = 3.5, y = 5) + 
  ylim(0, 5.5)

#test if chimp deletions between tiles predicts differential activity
chimp_loss_diffs <- data.frame(deletion_size = c(na.omit(all_active_df$deletion_size), na.omit(all_emvars_df$deletion_size)), 
                       Set = c(rep(x = "Active", length(na.omit(all_active_df$deletion_size))), rep(x = "Differentially Active", length(na.omit(all_emvars_df$deletion_size))))
) %>% dplyr::filter(deletion_size > 0)

wilcox.test(chimp_loss_diffs$deletion_size[which(chimp_loss_diffs$Set == "Active")], 
            chimp_loss_diffs$deletion_size[which(chimp_loss_diffs$Set == "Differentially Active")])

chimp_loss_activity <- ggplot(data = chimp_loss_diffs, aes(x = Set, y = deletion_size, fill = Set))+
  geom_boxplot(outliers = F) +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = "Size of chimp deletion(s) (bp)") + 
  scale_fill_manual(values = c("#56B4E9", "#0079B2")) + 
  annotate("text", label='NS', x = 1.5, y = 10) + 
  ylim(0, 11)

cor.test(all_emvars_df %>% dplyr::filter(align_type == "chimp_deletion", skew_sign_match) %>% pull(max_skew), 
         all_emvars_df %>% dplyr::filter(align_type == "chimp_deletion", skew_sign_match) %>% pull(deletion_size))

#plot skew against deletion size
deletion_size_plot <- ggplot(data = all_emvars_df %>% 
                               dplyr::filter(align_type == "chimp_deletion", skew_sign_match) %>% 
                               mutate(x_bins = cut(deletion_size, breaks = c(0,2,5,10,15,20,27))), 
                             aes(x = x_bins, y = max_skew)) + 
  geom_boxplot(fill = "lightgrey", outliers = F) +
  labs(x = "Size of chimp deletion(s) (base pairs)", y = "|Maximum log2 skew|") +
  scale_x_discrete(labels = c('1-2','3-5', '5-10','11-15', "15-20", "21-26")) + 
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) + 
  annotate("text", label='NS', x = 3.5, y = 5) + 
  ylim(0, 5.5)

#create the figure
chimp_loss_figure <- ggdraw() +
  draw_plot(chimp_loss_activity , x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(deletion_size_plot , x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(chimp_loss_activity_conserved, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(deletion_size_plot_conserved, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D"), 
                  size = 10,
                  x = c(0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.5, 0.5)) + 
  bgcolor("white")

ggsave(plot = chimp_loss_figure + panel_border(color = "black", size = 1), 
       filename = "~/Dropbox/Cartilage MPRA Paper/Code/fig_S7_chimp_losses.png", 
       device = "png", dpi = 300, height = 6.5, width = 6.5, units = "in")

