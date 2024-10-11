#!/usr/bin/env Rscript

################################################################################################################
## Filename: PopulationSummary.R
## Created: Oct 11, 2024
## Author(s): Muneeza Maqsood, Chung Kok
##
## Purpose: 
##      - to generate a population summary of amalgamated variants from all the samples in a batch
##          
## Instructions:  
##      - perform analysis on the SpliceHunter directory
##      - execute bash script 
##         
################################################################################################################

#Set working directory

setwd("./bams/tmp_bams/output/")


## load library
library(pacman)

p_load("tidyverse", 
       "Biostrings", 
       "BSgenome.Hsapiens.UCSC.hg19", 
       "dplyr", 
       "stringr",
       "randomForest",
       "DT", "ggplot2", "data.table")

pacman::p_load(readxl, tidyverse, xlsx, data.table, dplyr, purrr, matrixStats, ggplot2, openxlsx, magrittr, gapminder)
options(scipen = 999)


### Functions

#Filtering Hotspot genes from the HGNC_symbol column (33 genes):
hotspot_gene <- function(x){
  hotspotGene <-c( "ASXL1", "ATM","BCOR","BCORL1","BTG1","CALR","CDKN2A","CDKN2B","CUX1","DNMT3A","ETV6",
                   "EZH2","GATA2","IKZF1","JAK2","KDM6A","KIT","KMT2D","PAX5","PHF6","PPM1D","RB1",
                   "RUNX1","SETD1B","SETD2","SF3B1","SMC3","SRSF2","STAG2","TET2","TP53","UBE2A","WT1", "MTAP", "CDKN2A-AS1","CDKN2B-AS1", 
                   "EP300", "NUP98", "CREBBP", "RAD21", "RCBTB2")
  x <-  x %>% dplyr::filter(geneName %in% hotspotGene)
  x
}

###Filtering on read 
###To Create a matrix for comparison between more than 1 dfs ####

new_df <- function(x){
  
  cols_to_retain <- c("chr", "start", "stop", "startAnno", "stopAnno", "geneName", "del_size" , "intron_motif", "unique_count", "event",
                      "H_end", "H_start", "geneStrand", "distance_from_last5", "distance_from_first3", 
                      "strandDirection", "nExSkip", "exon", "diff_H_start_end", "delEx", "event5", "event3",
                      "impAnnot5", "impAnnot3", "JPKM")
  x <- x %>% dplyr::select(all_of(cols_to_retain)) %>%
    mutate(chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx_impAnnot5_impAnnot3 = paste(chr,start,
                                                                                                                                                          stop, startAnno, 
                                                                                                                                                          stopAnno, geneName,
                                                                                                                                                          del_size,intron_motif,
                                                                                                                                                          event, geneStrand,
                                                                                                                                                          strandDirection,nExSkip,
                                                                                                                                                          exon, delEx,
                                                                                                                                                          impAnnot5, impAnnot3,
                                                                                                                                                          sep = "__"),
           .before = chr)
  x
}

#Defining function to calculate stats (11) ## For comparing multiple datasets
calc_stat <- function(x){
  x1 <- x %>% remove_rownames() %>% 
    column_to_rownames(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx_impAnnot5_impAnnot3") %>%
    dplyr::select(starts_with("unique_count")) %>%
    type_convert(.) %>%
    transmute(`#PatientsWithVarCall` = rowSums(!is.na(.)))
  
  x2 <- x %>% remove_rownames() %>% 
    column_to_rownames(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx_impAnnot5_impAnnot3") %>%
    dplyr::select(starts_with("unique_count")) %>%
    type_convert(.) %>% rowwise() %>%
    transmute(median_count = round(median(c_across(everything()), na.rm =T), 0),
              mean_count = round(mean(c_across(everything()), na.rm =T), 2),
              SD_count = round(sd(c_across(everything()), na.rm = T), 2),
              CV_count = round((SD_count/mean_count), 1),
              max_count = max(c_across(everything(starts_with("unique"))), na.rm = T),
              min_count = min(c_across(everything(starts_with("unique"))), na.rm = T),
              total_count = sum(c_across(everything(starts_with("unique"))), na.rm = T),
              diff_max_min_count = abs(max_count - min_count),
              ratio_max_med_count = round((max_count/median_count), 0))
  
  x3 <- x %>% remove_rownames() %>% 
    column_to_rownames(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx_impAnnot5_impAnnot3") %>%
    dplyr::select(starts_with("JPKM")) %>%
    type_convert(.) %>% rowwise() %>%
    transmute(
      mean_JPKM = round(mean(c_across(everything()), na.rm =T), 2),
      SD_JPKM = round(sd(c_across(everything()), na.rm = T), 2),
      CV_JPKM = round((SD_JPKM/mean_JPKM), 1))
  
  
  x4 <- x %>% remove_rownames() %>% 
    column_to_rownames(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx_impAnnot5_impAnnot3") %>%
    dplyr::select(starts_with("H_end")) %>%
    type_convert(.) %>% rowwise() %>%
    transmute(max_H_end = max(c_across(everything(starts_with("H_end"))), na.rm = T))
  
  x5 <- x %>% remove_rownames() %>% 
    column_to_rownames(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx_impAnnot5_impAnnot3") %>%
    dplyr::select(starts_with("H_start")) %>%
    type_convert(.) %>% rowwise() %>%
    transmute(max_H_start = max(c_across(everything(starts_with("H_start"))), na.rm = T))
  
  x6 <- x %>% remove_rownames() %>%
    column_to_rownames(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx_impAnnot5_impAnnot3")
  x <- cbind(x6, x5, x4, x3, x2, x1)
}

del <- function(x){ #Calc deletion size
  x <- x %>% type_convert(.) %>%
    dplyr::mutate(del_size = abs(start - stop) )
  x
}

impAnno <- function(x){
  
  x <- x %>%
    dplyr::mutate(impAnnot5 = ifelse(abs(distance_from_last5) == 1, paste(event5, "start", sep = "-"), ifelse(abs(distance_from_last5) == 2, paste(event5, "end", sep = "-"), paste(event5, "middle", sep = "-") )),
                  impAnnot3 = ifelse(abs(distance_from_first3) == 1, paste(event3, "start", sep = "-"), ifelse(abs(distance_from_first3) == 2, paste(event3, "end", sep = "-"), paste(event3, "middle", sep = "-") ))) 
  x
  
}

#Removing Irrelevant cols (Not used by Sue in decision making)
rm_cols2 <- function(x){
  cols_to_remove <- c("mean_count", "SD_count", "diff_max_min_count", "SD_JPKM", "mean_JPKM")
  x <- x %>% 
    dplyr::select(-all_of(cols_to_remove))
  x
}


## Reading files:
all_samples <- list.files(pattern = "*_score.csv")

list_all_samples <- sapply(all_samples, read.csv, stringsAsFactors = F, simplify = F)

# filtering for gene of interest (33)
goi <- sapply(list_all_samples, hotspot_gene, simplify = F)

# intron motif annotation
list_all_samples2 <- sapply(goi, impAnno, simplify = F)

#Creating a df_list with equal number of columns:
shorterDF_list <- sapply(list_all_samples2, new_df, simplify = F)

#merging the list of dfs into a single df:
mergedDF <- bind_rows(shorterDF_list, .id = "id")

#pivot_wider according to the IDs of dataframes
mergedDF %>% 
  pivot_wider(id_cols = chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx_impAnnot5_impAnnot3,
              names_from = id,
              values_from = c(unique_count,H_end, H_start, diff_H_start_end, JPKM),
              names_vary = "slowest") -> df_wider

#Calculating some stats of the variants data and how many times are these variants been called (# of Patients & runs)
df_wider %>% calc_stat() %>% arrange(desc(`#PatientsWithVarCall`)) -> final_stat

# final_stat1 <- final_stat %>% calc_stat2()
final_stat1 <- final_stat %>% rm_cols2()
# Rearranging cols
final_stat2 <- final_stat1 %>% dplyr::select(`#PatientsWithVarCall`,  "max_count", "min_count", "median_count", "ratio_max_med_count", 
                                             "CV_count","CV_JPKM", "total_count", "max_H_start",
                                             "max_H_end", everything())
final_stat2 %>% rownames_to_column(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx_impAnnot5_impAnnot3") -> final_stat3
finale <- final_stat3 %>% separate(chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx_impAnnot5_impAnnot3, 
                                   c("chr","start", "stop", "startAnno", "stopAnno","geneName", "delSize", "intronMotif", "event", 
                                     "geneStrand", "strandDirection", "nExSkip", "exon", "delEx",
                                     "impAnnot5", "impAnnot3"),sep = "__") %>% suppressWarnings(type_convert(.))

fin3 <- finale %>% dplyr::select(-geneStrand, -strandDirection, -event, -total_count, -max_H_start, -max_H_end, -ratio_max_med_count)

#s Set working directory to extract the base directory name
setwd("../../../../.")

# Extract the folder name from the working directory
foldername <- basename(getwd())

# change working directory to where the final output file has to be stored
setwd("./SpliceHunter/bams/tmp_bams/output/NovelSplices/")

# Create a new workbook
wb <- createWorkbook()

# Add a worksheet to the workbook
addWorksheet(wb, "Sheet1")

# copy final df name to a new variable to write into the files
df <- fin3 %>% type_convert(.)

# Get the column names
column_names <- colnames(df)

# Define the color to apply to the columns starting with "unique"
color <- "#ccccff"  # You can change this to any desired color

# Iterate over the column names and add data to the worksheet
for (i in 1:length(column_names)) {
  column_name <- column_names[i]
  column_data <- df[[column_name]]
  
  # Write the column name and data to the worksheet
  hs <- createStyle(textDecoration = "BOLD", valign = "top", wrapText = T)
  writeData(wb, "Sheet1", column_name, startCol = i, startRow = 1, colNames = T, headerStyle = hs)
  writeData(wb, "Sheet1", column_data, startCol = i, startRow = 1, colNames = T, headerStyle = hs)
  
  if (grepl("^#Patients", column_name)) {
    cellStyle <- createStyle(fontColour = "blue")
    addStyle(wb, "Sheet1", cellStyle, rows = 1:nrow(df)+1, cols = i)
  }
  #Check if the column name starts with "unique" and apply the color
  if (grepl("^unique", column_name)) {
    cellStyle <- createStyle(fgFill =  color)
    addStyle(wb, "Sheet1", cellStyle, rows = 1:nrow(df)+1, cols = i)
  }
}

# Save the workbook to a file
writeDataTable(wb, "Sheet1", df, colNames = T, headerStyle = hs, rowNames = F)
saveWorkbook(wb, paste0(foldername, "_population_summary_splices.xlsx"), overwrite = TRUE)
