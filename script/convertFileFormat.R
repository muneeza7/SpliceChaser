#!/usr/bin/env Rscript

################################################################################################################
## Filename: convertFileFormat.R
## Created: October 11, 2024
## Author(s): Muneeza Maqsood, Chung Kok
##
## Purpose: 
##      - to convert file format of STAR junctions to bed format
##          
## Instructions:  
##      - make a new directory named "input" 
##      - will work on spliceJunctions_processed.tsv only
##         
################################################################################################################
rm(list=ls())

## load library
library(pacman)
p_load("tidyverse", "dplyr", "stringr")

## list all files from raw data
files <- list.files("./input/", pattern = "_processed.tsv", full.names = T)

## file convert to STAR junction output format
for(i in files) {
  id1 <- str_extract(i, "[a-zA-Z0-9-_\\. ]+_processed.tsv")
  fn <- str_trim(gsub(".spliceJunctions_processed.tsv|.star.b37", "", id1))

  print(paste("Processing...", fn, sep=""))
  dat <- read.table(i, header=T, sep=" ") %>%
    dplyr::select(chr, first_base_intron, last_base_intron, strand,
                  intron_motif, STAR_splice_type, unique_mappers, multi_mappers, 
                  max_splice_overhang) %>%
    mutate(chr=ifelse(str_detect(chr, "chr"),chr, paste0("chr",chr)))
  write.table(dat, file=paste0("./input/tmp/", fn, ".bed"), quote = F, sep = "\t", row.names = F,
              col.names = F)
}

print("Done...")
