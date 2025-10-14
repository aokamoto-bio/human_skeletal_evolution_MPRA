#create metadata tables for cartilage MPRA processing

#load functions from Tewhey lab
#download from their github and set up your own filepath #(https://github.com/tewhey-lab/MPRAmodel)
source("~/Desktop/Capellini_Lab/Weekly_Coding/MPRAmodel.R")

#load packages
library(Biostrings)

#load the required data

#read table as tab separated file
oligo_seqs <- read.table(file = "~/Desktop/Capellini_Lab/MPRA/Cartilage_MPRA/cartilage_MPRA_oligo_seqs2.tsv", header = T)

#clean up columns
cartilage_metadata <- oligo_seqs %>% 
  separate(col = name, sep = "_", into = c("region", "tile", "species"))
cartilage_metadata$ID <- paste(cartilage_metadata$region, cartilage_metadata$tile, sep ="_")
cartilage_metadata$species[which(is.na(cartilage_metadata$species))] <- "control"
cartilage_metadata$ID[which(cartilage_metadata$species == "control")] <- cartilage_metadata$region[which(cartilage_metadata$species == "control")]
cartilage_metadata$species <- gsub(pattern = "hs", replacement = "human", cartilage_metadata$species)

#add information about control type
cartilage_metadata$control_type <- NA
cartilage_metadata$control_type[which(cartilage_metadata$ID %in% MPRA_270bp_pos_controls$ID)] <- "Positive_Control"
cartilage_metadata$control_type[which(cartilage_metadata$ID %in% MPRA_270bp_neg_controls$ID)] <- "Negative_Control"

#add information about structural difference

#does the region map across multiple chimp chromosomes?
chimp_multi_match <- chimp_cartilage_regions_tiled_df %>% 
  dplyr::select(ID2) %>% 
  group_by(ID2) %>% 
  summarize(n = n()) %>% 
  dplyr::filter(n > 1) %>% 
  pull(ID2)

cartilage_metadata$chimp_multi_match <- "NO"
cartilage_metadata$chimp_multi_match[which(cartilage_metadata$ID %in% chimp_multi_match)] <- "YES"

#separate control data 
cartilage_metadata_controls <- cartilage_metadata %>% dplyr::filter(species == "control")
write_delim(cartilage_metadata_controls, file = "~/Desktop/Capellini_Lab/MPRA/Cartilage_MPRA/cartilage_metadata_controls.txt", delim = "\t", col_names = T)

#separate expr data 
cartilage_metadata_experimental <- cartilage_metadata %>% 
  dplyr::filter(species != "control") %>% 
  dplyr::group_by(ID) %>% 
  pivot_wider(names_from = species, values_from = sequence) %>% 
  separate(col = region, into = c("chr", "start"), sep = ":", remove = F)
cartilage_metadata_experimental$start <- as.numeric(cartilage_metadata_experimental$start)
cartilage_metadata_experimental$start[which(cartilage_metadata_experimental$tile == 2)] <- cartilage_metadata_experimental$start[which(cartilage_metadata_experimental$tile == 2)] + 270
cartilage_metadata_experimental$end <- (cartilage_metadata_experimental$start + 269)

#keep region information for generating beds
cartilage_metadata_experimental$region_start <- as.numeric(cartilage_metadata_experimental$start)
cartilage_metadata_experimental$region_start[which(cartilage_metadata_experimental$tile == 2)] <- cartilage_metadata_experimental$start[which(cartilage_metadata_experimental$tile == 2)] - 270

#find positions with specific overlaps
#return positions of dataset
#(input df should be able to generate a bed and include a 4th ID column)
overlap_MPRA_with_bed_set <- function(dataset, bed) {
  
  #create temporary bed file
  df_bed <- "~/Desktop/TEMP_ATAC_DIR/temp_MPRA_df.bed"
  write_bed(data = dataset, filename = df_bed, ncol = 4)
  
  #perform the overlap
  system(paste("bedtools intersect -a", df_bed, "-b", bed, "-wo > ~/Desktop/TEMP_ATAC_DIR/temp_MPRA_df_overlapped.bed"))
  
  #load overlapped df
  overlap_df <- read.delim(file = "~/Desktop/TEMP_ATAC_DIR/temp_MPRA_df_overlapped.bed", sep = "", head = F)
  
  #cleanup
  system(paste("rm", df_bed))
  system("rm ~/Desktop/TEMP_ATAC_DIR/temp_MPRA_df_overlapped.bed")
  
  #return the output
  return(overlap_df)
}

#overlap with HARs
HARs_overlap <- cartilage_metadata_experimental %>% 
  dplyr::select(chr, start, end, ID) %>% 
  overlap_MPRA_with_bed_set(bed = HARs)

#add HAR overlap data
cartilage_metadata_experimental <- HARs_overlap %>% 
  dplyr::select(V4, V8) %>% 
  rename(c(V4 = "ID", V8 = "HARs_Overlapped")) %>% 
  left_join(x = cartilage_metadata_experimental)

#overlap with hCONDELs
hCONDELs_overlap <- cartilage_metadata_experimental %>% 
  dplyr::select(chr, start, end, ID) %>% 
  overlap_MPRA_with_bed_set(bed = hCONDELs)

#add hCONDEL overlap data
cartilage_metadata_experimental <- hCONDELs_overlap %>% 
  dplyr::select(V4, V8) %>% 
  rename(c(V4 = "ID", V8 = "hCONDELs_Overlapped")) %>% 
  left_join(x = cartilage_metadata_experimental)

cartilage_metadata_experimental <- cartilage_metadata_experimental %>% 
  mutate(bp_diffs = compare.DNA(human, chimp))
cartilage_metadata_experimental$bp_diffs[which(cartilage_metadata_experimental$chimp_multi_match == "YES")] <- NA

colnames(cartilage_metadata_experimental)[5] <- "SNP"

HAQERs <- "/Users/alexanderokamoto/Desktop/Capellini_Lab/Structural_Variant_Beds/HAQERS_T2T_hg38.bed"
#overlap with HAQERs
HAQERs_overlap <- cartilage_metadata_experimental %>% 
  dplyr::select(chr, start, end, SNP) %>% 
  overlap_MPRA_with_bed_set(bed = HAQERs)

#add HAQERs overlap data
cartilage_metadata_experimental <- HAQERs_overlap %>% 
  dplyr::select(V4, V8) %>% 
  rename(c(V4 ="SNP", V8 = "HAQERs_Overlapped")) %>% 
  left_join(x = cartilage_metadata_experimental)

#fix with T2T data
cartilage_metadata_experimental$HAQERs_Overlapped <- NA
cartilage_metadata_experimental$HAQERs_Overlapped[which(cartilage_metadata_experimental$SNP %in% HAQERs_overlap$V4)] <- HAQERs_overlap$V8

#load dataframe with bed files path for ease
results_file <- "~/Desktop/Capellini_Lab/Weekly_Coding/Capellini_Beds/capellini_skeleton_heatmap_beds.txt"
results_df <- read.table(results_file, header = T)

#first, timepoint analysis
early_merged_bed <- "/Users/alexanderokamoto/Desktop/Capellini_Lab/Weekly_Coding/Capellini_Beds/early_merged_hg38.bed"
system(paste("cat", 
             paste(results_df$early_file_hg38[c(-29,-30)], collapse = " "), 
             "| cut -f'1-3' | sort -k 1,1 -k2,2n | bedtools merge -i stdin > ", 
             early_merged_bed, sep = " "))
late_merged_bed <- "/Users/alexanderokamoto/Desktop/Capellini_Lab/Weekly_Coding/Capellini_Beds/late_merged_hg38.bed"
system(paste("cat", 
             paste(results_df$late_file_hg38, collapse = " "), 
             "| cut -f'1-3' | sort -k 1,1 -k2,2n | bedtools merge -i stdin > ", 
             late_merged_bed, sep = " "))
all_merged_bed <- "/Users/alexanderokamoto/Desktop/Capellini_Lab/Weekly_Coding/Capellini_Beds/All_Capellini_merged_IDR_threshold0.05_hg38.bed"
early_only_bed <- "/Users/alexanderokamoto/Desktop/Capellini_Lab/Weekly_Coding/Capellini_Beds/early_ONLY_hg38.bed"
system(paste("bedtools subtract -a", all_merged_bed, "-b", late_merged_bed, "-A > ", early_only_bed, sep = " "))
late_only_bed <- "/Users/alexanderokamoto/Desktop/Capellini_Lab/Weekly_Coding/Capellini_Beds/late_ONLY_hg38.bed"
system(paste("bedtools subtract -a", all_merged_bed, "-b", early_merged_bed, "-A > ", late_only_bed, sep = " "))
tp_shared_bed <- "/Users/alexanderokamoto/Desktop/Capellini_Lab/Weekly_Coding/Capellini_Beds/tp_shared_hg38.bed"
system(paste("bedtools subtract -a", all_merged_bed, "-b", early_only_bed, "-A | bedtools subtract -a stdin -b", late_only_bed , " -A > ", tp_shared_bed, sep = " "))

#add timepoint info to metadata
cartilage_metadata_experimental$timepoint <- "NA"
early_regions <- cartilage_metadata_experimental %>% ungroup() %>% 
  dplyr::select(chr, start, end, region) %>% 
  overlap_MPRA_with_bed_set(bed = early_only_bed) %>% 
  dplyr::pull(V4) %>% unique()
cartilage_metadata_experimental$timepoint[which(cartilage_metadata_experimental$region %in% early_regions)] <-"EARLY_ONLY"
late_regions <- cartilage_metadata_experimental %>% 
  ungroup() %>% 
  dplyr::select(chr, start, end, region) %>% 
  overlap_MPRA_with_bed_set(bed = late_only_bed) %>% 
  dplyr::pull(V4) %>% 
  unique()
cartilage_metadata_experimental$timepoint[which(cartilage_metadata_experimental$region %in% late_regions)] <-"LATE_ONLY"

shared_regions <- cartilage_metadata_experimental %>% 
  ungroup() %>% 
  dplyr::select(chr, start, end, region) %>% 
  overlap_MPRA_with_bed_set(bed = tp_shared_bed) %>% 
  dplyr::pull(V4) %>% 
  unique()
cartilage_metadata_experimental$timepoint[which(cartilage_metadata_experimental$region %in% shared_regions)] <-"SHARED"

#now the tissue overlaps
results_df$tp_merged_hg38 <- gsub(x = results_df$unique_hg38, 
                                  pattern = "_BOTH_TP_UNIQUE_hg38.bed", 
                                  replacement = "_BOTH_TP_hg38.bed")
results_df$tp_merged_hg38[29:30] <- results_df$late_file_hg38[29:30]

for(x in 1:nrow(results_df)){
  if(!x %in% 29:30){
    system(paste("cat ", results_df$early_file_hg38[x], " ", results_df$late_file_hg38[x], " | cut -f'1-3' | sort -k 1,1 -k2,2n | bedtools merge -i stdin > ", results_df$tp_merged_hg38[x], sep ="")) 
 
  }
}


MPRA_multiinter <- correct_multiinter(results_df$tp_merged_hg38, save_bed = T, output_bed_name = "~/Desktop/290k_MPRA_Cluster_Results/results/MPRA_multiinter_tissue_count.bed")
colnames(MPRA_multiinter) <- c("chr", "start", "end", "num", "list", results_df$tissue)

cartilage_metadata_experimental <- cartilage_metadata_experimental %>% ungroup()

cartilage_metadata_experimental %>% 
  dplyr::select(chr, start, end, region) %>% 
  write_bed(filename = "~/Desktop/290k_MPRA_Cluster_Results/results/cartilage_MPRA_metadata_hg38.bed", ncol = 4)

#perform the overlap
system(paste("bedtools intersect -a", 
"~/Desktop/290k_MPRA_Cluster_Results/results/cartilage_MPRA_metadata_hg38.bed", 
"-b", 
"~/Desktop/290k_MPRA_Cluster_Results/results/MPRA_multiinter_tissue_count.bed",
"-wo > ~/Desktop/290k_MPRA_Cluster_Results/results/metadata_MPRA_multiinter_tissue_count.bed"))

#load overlapped df
MPRA_overlap_df <- read.delim(file = "~/Desktop/290k_MPRA_Cluster_Results/results/metadata_MPRA_multiinter_tissue_count.bed", sep = "", head = F)
colnames(MPRA_overlap_df) <- c("chr", "start", "end", "region", "chr.ATAC", "start.ATAC", "end.ATAC", "tissue_overlap_n", "list", results_df$tissue, "Overlap_len")
MPRA_overlap_df2 <- MPRA_overlap_df %>% 
  dplyr::select(-end, -chr.ATAC, -start.ATAC, -end.ATAC, -Overlap_len)

#some sequences overlapped two regions
problemmatic_regions <- c("chr1:103078499_1", "chr14:105477429_1", "chr14:105477429_2", "chr2:33068197_2", "chr5:94569545_1", "chrX:55183023_2", "chrX:55183259_1")

#fix the length match question
cartilage_metadata_reduced <- cartilage_metadata_experimental[, 1:14]
cartilage_metadata_reduced$align_type <- NA
cartilage_metadata_reduced$align_type[which(cartilage_metadata_reduced$SNP %in% chimp_cartilage_regions_tiled_df$ID2[which(chimp_cartilage_regions_tiled_df$width == 270)])] <- "SNPs_only"
chimp_cartilage_regions_tiled_df$ID2[which(chimp_cartilage_regions_tiled_df$width ==270)]

#now for chimp gap but same length
cartilage_metadata_reduced$align_type[which(cartilage_metadata_reduced$SNP %in% 
                                              setdiff(cartilage_regions_tiled_df_filtered$ID2[which(cartilage_regions_tiled_df_filtered$width_chimp == 270)], cartilage_metadata_reduced$SNP[which(cartilage_metadata_reduced$align_type == "SNPs_only")]))] <- "chimp_gap_simple"

cartilage_metadata_reduced$bp_diffs[which(is.na(cartilage_metadata_reduced$align_type))] <- NA
cartilage_metadata_reduced$bp_diffs[which(cartilage_metadata_reduced$align_type == "chimp_gap_simple")] <- cartilage_metadata_reduced %>% 
  dplyr::filter(align_type == "chimp_gap_simple") %>% 
  rowwise() %>% 
  mutate(bp_diffs = compare.DNA(human, chimp)) %>% 
  pull(bp_diffs)

cartilage_metadata_reduced$align_type[which(cartilage_metadata_reduced$SNP %in% 
                                              cartilage_regions_tiled_df_filtered2$ID2[which(cartilage_regions_tiled_df_filtered2$padding > 0)])] <- "chimp_deletion"

cartilage_metadata_reduced <- cartilage_regions_tiled_df_filtered2 %>% 
  mutate(SNP = ID2, deletion_size = padding) %>% 
  dplyr::select(SNP, deletion_size) %>% 
  merge(x= cartilage_metadata_reduced, all.y = F)
cartilage_metadata_reduced$deletion_size <- NULL

cartilage_metadata_reduced <- cartilage_metadata_reduced %>% 
  mutate(score = pairwiseAlignment(human, chimp, scoreOnly =T))
length(which(cartilage_metadata_reduced$score > 0))
length(which(cartilage_metadata_reduced$score <= 0))

misaligned_SNPs <- cartilage_metadata_reduced$SNP[which(cartilage_metadata_reduced$score <= 0)]

#remove misaligned SNPs
cartilage_metadata_experimental <- cartilage_metadata_experimental %>% 
  dplyr::filter(!SNP %in% misaligned_SNPs)

cartilage_metadata_experimental$bp_diffs <- NULL
cartilage_metadata_experimental <- cartilage_metadata_reduced %>% 
  dplyr::select(SNP, bp_diffs, align_type, deletion_size, score) %>% 
  merge(x = cartilage_metadata_experimental, by = "SNP")

#saving and load the current version
write_delim(cartilage_metadata_experimental, file = "~/Desktop/Capellini_Lab/MPRA/Cartilage_MPRA/cartilage_metadata_experimental_alignable.txt", delim = "\t", col_names = T)
cartilage_metadata_experimental <- read.table(file = "~/Desktop/Capellini_Lab/MPRA/Cartilage_MPRA/cartilage_metadata_experimental_alignable.txt", header = T)
