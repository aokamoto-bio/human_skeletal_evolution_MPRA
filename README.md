

**Code for the preprint:** 

Alexander Okamoto,  Clarissa Coveney, Danalaxshmi S Ganapathee, & Terence Capellini. Massively parallel functional screen identifies thousands of regulatory differences in human versus chimpanzee postcranial skeletal development

DOI COMING SOON

**Note**

Initial processing of the sequencing data was performed on the Harvard FAS computing cluster. Other parts of this script were written for a smaller, internal lab computing cluster. It may not work on a personal computer.

The code provided here is exactly what was run to perform these analyses but may not run on another system without modifications to ensure proper software, computational power, and directory structures. Please reach out to me at aokamoto@g.harvard.edu if you wish to replicate these analyses and the provided code is insufficient. I will try to help you to the best of my abilities.

Some functions used in the study are available from: https://github.com/aokamoto-bio/Human_Autopod_Evolution/blob/main/Custom_Analysis_Functions.R

**Files:**

1. Cartilage_MPRA_Design_270bp.Rmd details the construction of the library based on ATAC-seq data. 
2. Cartilage_MPRA_Processing.Rmd contains the commands used to run the MPRAsuite tools MPRAmatch and MPRA count to map oligos to barcodes and then count the barcodes.
3. 290k_Cartilage_Build_Attributes.R create attributes file for MPRAmodel
4. 290k_Cartilage_MPRA_Processing_Short.R run MPRAmodel to determine active and differentially active tiles
5. 290k_Cartilage_MPRA_Metadata.R code for constructing metadata dataframe for each oligo pair
6. MPRA_1_Controls.R for analysis of control oligos
7. MPRA_2_Activity_Profiling.R for analysis of active tiles
8. MPRA_3_Differential_Activity_Profiling.R for analysis of differentially active tiles
9. MPRA_4_Evolutionary_Testing.R for analysis of evolutionary signals in the data
10. MPRA_5_chimp_losses.R containing analysis of chimp losses for Figure S7
