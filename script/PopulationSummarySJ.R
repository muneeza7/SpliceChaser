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
##      - perform analysis on the SpliceChaser directory
##      - execute bash script 
##         
################################################################################################################

#Set working directory

setwd("./bams/tmp_bams/output/")


## load library
library(pacman)

pacman::p_load(tidyverse, data.table, dplyr, openxlsx)
options(scipen = 999)


### Functions

###Filtering on H_end & H_start
div_read_fil <- function(x){
  x <- x %>% dplyr::filter(H_end >= 2.67 & H_start >=2.67)
  x$JRPM <- round(x$JRPM, 0)
  x
###To Create a matrix for comparison between more than 1 dfs ####

new_df <- function(x){
  
  cols_to_retain <- c("chr", "start", "stop", "startAnno", "stopAnno", "geneName", "del_size" , "intron_motif", "unique_count", "event",
                      "H_end", "H_start", "geneStrand", "distance_from_last5", "distance_from_first3", 
                      "strandDirection", "nExSkip", "exon", "diff_H_start_end", "delEx", "event5", "event3", "JRPM")
  x <- x %>% dplyr::select(all_of(cols_to_retain)) %>% 
    mutate(chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx = paste(chr,start,
                                                                                                                                                          stop, startAnno, 
                                                                                                                                                          stopAnno, geneName,
                                                                                                                                                          del_size,intron_motif,
                                                                                                                                                          event, geneStrand,
                                                                                                                                                          strandDirection,nExSkip,
                                                                                                                                                          exon, delEx,
                                                                                                                                                          sep = "__"),
           .before = chr)
  x
}


#Defining function to calculate stats (11) ## For comparing multiple datasets
calc_stat <- function(x){
  x1 <- x %>% remove_rownames() %>% 
    column_to_rownames(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx") %>%
    dplyr::select(starts_with("unique_count")) %>%
    type_convert(.) %>%
    transmute(`#PatientsWithVarCall` = rowSums(!is.na(.)))
  
  x2 <- x %>% remove_rownames() %>% 
    column_to_rownames(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx") %>%
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
    column_to_rownames(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx") %>%
    dplyr::select(starts_with("JRPM")) %>%
    type_convert(.) %>% rowwise() %>%
    transmute(
      mean_JRPM = round(mean(c_across(everything()), na.rm =T), 2),
      SD_JRPM = round(sd(c_across(everything()), na.rm = T), 2),
      CV_JRPM = round((SD_JRPM/mean_JRPM), 1))
  
  
  x4 <- x %>% remove_rownames() %>% 
    column_to_rownames(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx") %>%
    dplyr::select(starts_with("H_end")) %>%
    type_convert(.) %>% rowwise() %>%
    transmute(max_H_end = max(c_across(everything(starts_with("H_end"))), na.rm = T))
  
  x5 <- x %>% remove_rownames() %>% 
    column_to_rownames(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx") %>%
    dplyr::select(starts_with("H_start")) %>%
    type_convert(.) %>% rowwise() %>%
    transmute(max_H_start = max(c_across(everything(starts_with("H_start"))), na.rm = T))
  
  x6 <- x %>% remove_rownames() %>%
    column_to_rownames(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx")
  x <- cbind(x6, x5, x4, x3, x2, x1)
}

#Calc deletion size
del <- function(x){
  x <- x %>% type_convert(.) %>%
    dplyr::mutate(del_size = abs(start - stop) )
  x
}

#Removing Irrelevant cols 
rm_cols2 <- function(x){
  cols_to_remove <- c("mean_count", "SD_count", "diff_max_min_count", "SD_JPKM", "mean_JPKM")
  x <- x %>% 
    dplyr::select(-all_of(cols_to_remove))
  x
}


## Reading files:
all_samples <- list.files(pattern = "*_score.csv")

list_all_samples <- sapply(all_samples, read.csv, stringsAsFactors = F, simplify = F)

# H_score filtering
list_all_samples1 <- sapply(list_all_samples, div_read_fil, simplify = F)

#Creating a df_list with equal number of columns:
shorterDF_list <- sapply(list_all_samples1, new_df, simplify = F)

#merging the list of dfs into a single df:
mergedDF <- bind_rows(shorterDF_list, .id = "id")

#pivot_wider according to the IDs of dataframes
mergedDF %>% 
  pivot_wider(id_cols = chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx,
              names_from = id,
              values_from = c(unique_count,H_end, H_start, diff_H_start_end, JRPM),
              names_vary = "slowest") -> df_wider

#Calculating some stats of the variants data and how many times are these variants been called (# of Patients & runs)
df_wider %>% calc_stat() %>% arrange(desc(`#PatientsWithVarCall`)) -> final_stat

# final_stat1 <- final_stat %>% calc_stat2()
final_stat1 <- final_stat %>% rm_cols2()
# Rearranging cols
final_stat2 <- final_stat1 %>% dplyr::select(`#PatientsWithVarCall`,  "max_count", "min_count", "median_count", "ratio_max_med_count", 
                                             "CV_count","CV_JRPM", "total_count", "max_H_start",
                                             "max_H_end", everything())
final_stat2 %>% rownames_to_column(var = "chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx") -> final_stat3
finale <- final_stat3 %>% separate(chr_start_stop_startAnno_stopAnno_geneName_delSize_intronMotif_event_geneStrand_strandDirection_nExSkip_exon_delEx, 
                                   c("chr","start", "stop", "startAnno", "stopAnno","geneName", "delSize", "intronMotif", "event", 
                                     "geneStrand", "strandDirection", "nExSkip", "exon", "delEx"),sep = "__") %>% suppressWarnings(type_convert(.))

fin3 <- finale %>% dplyr::select(-geneStrand, -strandDirection, -event, -total_count, -max_H_start, -max_H_end, -ratio_max_med_count)

#Set working directory to extract the base directory name
setwd("../../../../.")

# Extract the folder name from the working directory
foldername <- basename(getwd())

# change working directory to where the final output file has to be stored
setwd("./SpliceChaser/bams/tmp_bams/output/NovelSplices/")

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
