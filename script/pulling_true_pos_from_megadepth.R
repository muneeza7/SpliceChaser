#!/usr/bin/env Rscript

################################################################################################################
## Filename: pulling_true_pos_from_megadepth.R
## Created: Oct 11, 2024
## Author(s): Muneeza Maqsood, Chung Kok
##
## Purpose: 
##      - to pull extra variants drectly from bams that are missed by STAR due to mapping parameters
##          
## Instructions:  
##      - perform analysis on the SpliceHunter directory
##      - put all input files (SJ.out.tab) into input folder 
##      - execute bash script 
##         
################################################################################################################
rm(list=ls())


## load library
library(pacman)

p_load("tidyverse", 
       # "Biostrings", 
       # "BSgenome.Hsapiens.UCSC.hg19", 
       "dplyr", 
       "stringr")

p_load("tidyverse", 
       "Biostrings",
       "BSgenome.Hsapiens.UCSC.hg19",
       "dplyr", 
       "stringr",
       "annotate",
       "org.Hs.eg.db",
       "TxDb.Hsapiens.UCSC.hg19.knownGene",
       "GenomicRanges",
       "Homo.sapiens",
       "GenomicFeatures",
       "GenomicAlignments",
       "SummarizedExperiment",
       "MatrixGenerics")

## list files
md <- list.files(path="./bams/tmp_bams/", pattern="_jxs.tsv", full.names = T)
star <- list.files(path="./input/", pattern="SJ.out.tab", full.names = T)

for(k in 1:length(md)) {
  
  # reading in file
  dat <- read.table(md[k], sep="\t", header=FALSE, stringsAsFactors = FALSE)
  
  ## extract id
  id1 <- str_extract(md[k], "[a-zA-Z0-9-_\\.]+.bam.all_jxs.tsv")
  id <- gsub(".bam.all_jxs.tsv|X__", "", id1)
  
  print(paste("processing ...  ", id1, sep=""))
  
  # Filtering reads that have deleted sequence or N in the CIGAR string
  dat <- dat %>% 
    filter(V7=="1") %>%
    mutate(n_count = str_count(V6, "N"),
           m_count=str_count(V6, "M")) %>%
    filter(n_count<=2 & m_count<=3) %>%
    mutate(n_skipped=(abs(V3-V4)+1))
  
  # Filtering the reads with deletions ranging between 15-60 bp
  dat <- dat %>%
    group_by(V2, V3, V4, V5, V7, n_skipped) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    mutate(n_skipped=as.numeric(n_skipped)) %>%
    filter(n>15 & (n_skipped >=15 & n_skipped<=60))
  
  # storing the remaining variants in a dataframe in the same format as STAR's SJ.out.tab file
  dat3 <- dat %>%
    mutate(n_skipped=factor(n_skipped)) %>%
    group_by(V2,V3,V4, n_skipped) %>%
    summarise(n1=n(), total_reads=sum(n)) %>%
    ungroup() %>%
    mutate(n_skipped=as.numeric(as.character(n_skipped))) %>%
    dplyr::select(-n1) %>%
    mutate(row_id = paste(V2, V3, V4, sep = "_"), .before = V2) %>%
    as.data.frame()
  
  dat4 <- dat3 %>% transmute(row_id = row_id, 
                             V1 = dat3$V2,
                             V2 = dat3$V3,
                             V3 = dat3$V4,
                             V4 = "0",
                             V5 = "0",
                             V6 = "0",
                             V7 = total_reads,
                             V8 = "0",
                             V9 = n_skipped
  )
  # Reading STAR's SJ.out.tab files to check if the shortlisted variants from above are reported by STAR or not
  star_splice <- read.table(star[k], sep="\t", header=FALSE, stringsAsFactors = FALSE)
  star_splice %>% mutate(row_id = paste(V1, V2, V3, sep = "_"), .before = V1) -> ss
  
  # Getting the variants that were not called by STAR but were detected by Megadepth in the above processing
  as.data.frame(setdiff(dat4$row_id, ss$row_id)) -> tst
  
  # Appending STAR's output splice junction file with newly detected variants from Megadepth
  if(length(rownames(tst)) == 0){
    comb <- ss
  }else{
    for (r in 1:length(rownames(tst))){
      tst[r,] -> hit
      which(dat4$row_id == hit) -> pos
      dat5 <- dat4[pos, ]
      rbind(ss, dat5) -> ss
      ss -> comb
    }
  }
  
  comb1 <- comb %>%
    dplyr::rename(
      chr = V1,
      first_base_intron = V2,
      last_base_intron = V3,
      strand = V4,
      intron_motif = V5,
      STAR_splice_type = V6,
      unique_mappers = V7,
      multi_mappers = V8,
      max_splice_overhang = V9) %>%
    dplyr::select(-row_id) %>%
    dplyr::mutate(strand = ifelse(strand==1, paste("+"), ifelse(strand==2, paste("-"), paste("0"))))
  
  # Writing into file
  write.table(comb1, file = paste("./input/",
                                  id, "_spliceJunctions_processed.tsv"),row.names = F)
}
