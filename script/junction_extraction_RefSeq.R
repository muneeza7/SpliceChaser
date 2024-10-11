#!/usr/bin/env Rscript

################################################################################################################
## Filename: junction_extraction_RefSeq.R
## Created: October 11, 2024
## Author(s): Muneeza Maqsood, Chung Kok
##
## Purpose: 
##      - to identify cryptic and atypical splicing
##      - to detect splice inferred inter- and intra-genic deletions
##          
## Instructions:  
##      - perform analysis on the SpliceHunter directory
##      - execute bash script 
##         
################################################################################################################
rm(list=ls())


## load library
library(pacman)

p_load("stringr","tidyverse", 
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


# Make the annotating function. It will annotate the intervals with gene_Ids
annotateIntervals <-function(intervals, txdb){
  stopifnot(is(intervals, "GRanges"), is(txdb, "TxDb"))
  anno = genes(txdb)
  olaps = findOverlaps(intervals, anno)
  mcols(olaps)$gene_id = genes$gene_id[subjectHits(olaps)]
  intervals_factor = factor(queryHits(olaps), levels=seq_len(queryLength(olaps)))
  intervals$gene_id = splitAsList(mcols(olaps)$gene_id, intervals_factor)
  intervals
}

# Using the reference fasta and fai files (we use for mapping) I created a .bed file with Chr info and lengths and removed other
# columns except the 1st 3 and saved the file as .csv
### Creating txdb object to annotate the coordinates:
# Reading the bed file:
bed <- read.table("./database/WholeGenomeFasta/human_hg19.csv", header = F, sep = ",")

# Used the bed file info to generate chrominfo obj
chrominfo <- data.frame(chrom = bed[,1],
                        length = bed[, 3]) %>% mutate(is_circular = FALSE)
# Copied RefSeq_hg19 url from Google (NCBI)
metadata <- data.frame(name = "Resource URL",
                       value = paste0("https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.25/"))
# Putting the objects created above in the following function to create the txdb obj
txdb <- makeTxDbFromGFF(file = "./database/RefSeq.gtf",
                        format = "gtf",
                        chrominfo = chrominfo,
                        dataSource = "entrez",
                        organism = "Homo sapiens",
                        metadata = metadata)

# Make an header for the data.frame "myIntervals" containing
# the coordinates files
header <- c("Chromosome", "Start", "End")

# column 1: chromosome
# column 2: first base of the intron (1-based)
# column 3: last base of the intron (1-based)
# column 4: strand (0: undefined, 1: +, 2: -)
# column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:
#   AT/AC, 6: GT/AT
# column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
# column 7: number of uniquely mapping reads crossing the junction
# column 8: number of multi-mapping reads crossing the junction
# column 9: maximum spliced alignment overhang

## list files
files <- list.files(path="./input/tmp/", pattern="_mapped.bed", full.names = T)

## load reference
load("./database/GTEx ref v7 junction database.RData")

# load reference without intron for transcript/gene length
geneLen <- read.table("./database/RefSeq_ref_without_Intron_v1.bed")
geneLen <- geneLen %>% dplyr::select(V5, V8) %>%
  dplyr::rename(t_len = V8) %>% unique()
## load gene_count files
gf <- list.files(path = "./GeneCounts/", pattern = "ReadCounts.tsv", full.names = T)

## gene of interest
goi <- c( "ASXL1", "ATM", "BCOR", "BCORL1", "BTG1", "CALR",
          "CBL", "CDKN2A", "CDKN2B", "CEBPA", "CREBBP", "CRLF2", "CSF3R",
          "CUX1", "DNMT3A", "ETV6", "EED", "EP300", "EZH2", "FAT1", "FAT3",
          "FBXW7", "FLT3", "GATA2", "IKZF1", "JAK2",
          "KDM5A", "KDM6A", "KIT", "KMT2D", "MPL", "NCOR1", "NF1",
          "NOTCH1", "NPM1", "NSD2", "NT5C2", "PAX5", "PHF6", "PPM1D",
          "PTEN", "PTPN11", "RAD21", "RB1", "RUNX1", "SETBP1", "SETD1B",
          "SETD2", "SF3B1", "SH2B3", "SMARCA4", "SMC1A", "SMC3", "SRSF2",
          "STAG2", "STAT5B", "TET2", "TP53", "UBE2A", "WHSC1",
          "WT1", "ZRSR2", "MTAP", "EP300", "NUP98", "CREBBP", "RAD21")

for(k in 1:length(files)) {
  
  # reading in file
  dat <- read.table(files[k], sep="\t", header=FALSE, stringsAsFactors = FALSE)
  
  ## extract id
  id1 <- str_extract(files[k], "[a-zA-Z0-9-_\\.]+_mapped.bed")
  id <- gsub("_mapped.bed|X__", "", id1)
  
  print(paste("processing ...  ", id1, sep=""))
  
  ## total splicing events unique counts
  totalDat <- dat %>%
    distinct(V1, V2, V3, .keep_all = T) %>%
    summarise(totalCounts=sum(V7))
  
  ## Reading Star's output gene counts
  gdat <- read.table(gf[k], sep = "\t", header = F, stringsAsFactors = F)
  gdat <- gdat[-c(1:4),] %>% remove_rownames() %>%
    dplyr::select(-V2, -V4)
  
  ## library size/total read counts
  totalGdat <- gdat %>% summarise(totalGCounts=sum(V3))
  
  ## Calculating "per million" factor
  pm <- as.numeric(totalGdat)/1000000
  
  # ## RPM Calc
  gdat2 <- gdat %>% mutate(rpm=V3/as.numeric(pm)) %>%
    dplyr::rename(Symbol=V1, Scounts=V3)
  
  
  geneDat2 <- dat %>%
    group_by(V18) %>%
    summarise(counts=sum(V7)) %>%
    dplyr::rename(Symbol=V18) 
  
  ## SpliceHunter
  posdat.1 <- dat %>%
    mutate(V4.1=V4) %>%
    mutate(V4=ifelse(V4=="+" | V4=="-", V4, V13)) %>%
    arrange(V1, V2, V3, V4, V18, V11, V12, V17) %>%
    mutate(overlap=V19-V20) %>%
    dplyr::filter(! (V15=="intron" & overlap==0)) %>% 
    group_by(V1, V2, V3, V14) %>%
    dplyr::mutate(
      prop_uniq_multi=V7/(V7+V8),
      first5 = ifelse(V13=="+", dplyr::first(V11), dplyr::last(V12) ), ## include strand ifelse
      last5 = ifelse(V13=="+", dplyr::first(V12), dplyr::last(V11)),
      first3 = ifelse(V13=="+", dplyr::last(V11), dplyr::first(V12)),
      last3 = ifelse(V13=="+", dplyr::last(V12), dplyr::first(V11)),
      event5 = ifelse(V13=="+", dplyr::first(V15), dplyr::last(V15)),
      event3 = ifelse(V13=="+", dplyr::last(V15), dplyr::first(V15)),
      location5 = ifelse(V13=="+", dplyr::first(V17), dplyr::last(V17)),
      location3 = ifelse(V13=="+", dplyr::last(V17), dplyr::first(V17)),
      strandDirection= ifelse(V13==V4, "same", "opposite")) %>%
    mutate(event=paste(event5, event3, sep="-"),
           nExSkip=ifelse(abs(location5-location3)==0, abs(location5-location3), abs(abs(location5-location3)-1)),
           exon=paste(event5, location5, event3, location3, sep="")) %>% ## collapse exon number per transcript
    mutate(exon=gsub("exon", "E", exon),
           exon=gsub("intron", "I", exon),
           withinEx5=ifelse(V13=="+" & V2>first5 & V2 < last5, "Yes", 
                            ifelse(V13=="-" & V3<first5 & V3>last5, "Yes", "No")),
           withinEx3=ifelse(V13=="+" & V3>first3 & V3 < last3, "Yes", 
                            ifelse(V13=="-" & V2<first3 & V2>last3, "Yes", "No"))) %>%
    ungroup() %>%
    group_by(V1, V2, V3) %>%
    mutate(nGene = length(unique(V18)),
           nTranscripts=length(unique(V14)), 
           distance_from_last5= ifelse(V4=="+", V2-last5, V3-last5),
           distance_from_first3= ifelse(V4=="+", V3-first3, V2-first3)) %>%
    mutate(total=paste(event, collapse = "|"),
           Gene=paste(unique(V18), collapse="|")) %>%
    mutate(abnormality=paste(ifelse(V4=="+" & str_detect(total, "exon-exon") & abs(distance_from_last5)==1 & 
                                      abs(distance_from_first3)==2 & nExSkip<1 & nGene==1 & strandDirection=="same", "normal",
                                    ifelse(V4=="-" & str_detect(total, "exon-exon") & 
                                             abs(distance_from_last5)==2 & abs(distance_from_first3)==1 & 
                                             nExSkip<1 & nGene==1 & strandDirection=="same", "normal", "abnormal")), 
                             collapse="|")) %>%
    as.data.frame()
  
  ### filter from GTEx reference
  ref1 <- ref %>% 
    dplyr::select(name, fraction)

  
  posdat1 <- posdat.1 %>% 
    mutate(name=paste(V1,V2,V3, sep="_")) %>%
    dplyr::filter(V7 > 0) %>%
    left_join(ref1, by=c("name"="name")) %>%
    dplyr::rename(GTExFraction=fraction) %>%
    dplyr::select(-name)
  
  
  ## filter no of reads and choose largest total exon num transcript for reporting
  bb4 <- posdat1 %>%
    group_by(V1, V2, V3, V18) %>%
    slice_max(V16, with_ties = F) %>%
    dplyr::filter((str_detect(abnormality, "\\babnormal\\b")) | GTExFraction < 0.01) %>%
    left_join(geneDat2, by=c("V18"="Symbol")) %>%
    mutate(JunctionRatio = V7/counts,
           JRPM = (V7 * 10^6)/as.integer(counts)) %>%
    left_join(geneLen, by=c("V14"="V5")) %>%
    mutate(JPKM=round(JRPM/(t_len/1000), 3)) %>%
    mutate(delEx=ifelse(withinEx5=="Yes" & withinEx3=="Yes" & nExSkip == "0", paste0("del_", abs(V3-V2+1), "bp"), ## need to include exon-exon or !intron-intron
                        ifelse(withinEx5=="Yes" & withinEx3=="Yes" & nExSkip >= 1, paste0("del", location5 +1,"-", location3 -1),
                               ifelse(withinEx5=="No" & withinEx3=="Yes" & nExSkip >= 1, paste0("del", location5,"-", location3-1),
                                      ifelse(withinEx5=="Yes" & withinEx3=="No" & nExSkip >= 1, paste0("del", location5+1,"-", location3),
                                             ifelse(withinEx5=="No" & withinEx3=="No" & nExSkip >= 1, paste0("del", location5,"-", location3),
                                                    ifelse(nGene > 1, "del_>1genes",
                                                           ifelse(!str_detect(event, "exon-exon"), "intron",
                                                                  ifelse(abs(distance_from_last5) >2 | abs(distance_from_first3) >2, "non_exon_boundary", ""))))))))) %>%
    
    mutate(normalJunction=str_extract(abnormality, "\\bnormal\\b")) %>%
    mutate(backSpliceInfo=ifelse(V4!=V13 & nGene==1, "backSplice", NA)) %>%
    group_by(V1, V2, V3, .keep_all = T) %>%
    as.data.frame()
 
  ##junction burden per sample
  
  # Junctions considered for the 'junction burden' calculation were all tumor sample 
  # junctions not found in core normal samples. The total junction count per patient 
  # was divided by the mapped read count of the sample divided by 10 000 
  # (scaling to 'per Mb' with the assumption of 100-bp reads) to
  # give the final junction burden. 
  
  totalJB <- bb4 %>%
    distinct(V1, V2, V3, .keep_all = T) %>%
    summarise(totalJ=sum(V7)) %>%
    mutate(JBpS= totalJ * 10000/(as.integer(totalDat)),
           sampleName=id,
           total_mapped=totalDat) %>%
    dplyr::select(sampleName, totalJ, total_mapped, JBpS)
  write.csv(totalJB, file=paste("./results/", id, "_junction burden per sample.csv", sep=""), row.names = F)
  
  ## final cleanning up
  final <- bb4 %>% 
    dplyr::filter( (V18 %in% goi) ) %>%
    dplyr::rename(chr=V1,
                  start=V2,
                  stop=V3,
                  readStrand=V4.1,
                  intron_motif=V5,
                  unique_count = V7,
                  multimapped_count=V8,
                  max_overhang = V9,
                  geneStrand=V13,
                  transcriptName=V14,
                  type=V15,
                  totalEx=V16,
                  ExNum=V17,
                  geneName=V18,
                  V42 = .keep_all) %>%
    mutate(intron_motif=case_when(
      intron_motif == 0 ~ "non-canonical",
      intron_motif == 1 ~ "GT/AG",
      intron_motif == 2 ~ "CT/AC",
      intron_motif == 3 ~ "GC/AG",
      intron_motif == 4 ~ "CT/GC",
      intron_motif == 5 ~ "AT/AC",
      intron_motif == 6 ~ "GT/AT")) %>% 
    dplyr::select(-V6, -V10, -V11, -V12, -V19, -V20, -V4, -V42,
                  -overlap, -total, -abnormality, -ExNum, -type) %>%
    mutate(sampleName=id) %>%
    as.data.frame()
  
  ## In case the sample have all the normal variants and nothhing passed after filtering
  if (nrow(final)==0){
    write.csv(paste("WARNING - NO RESULTS FOUND!"), file=paste("./results/", id,"_result_final.csv", sep=""), row.names = F)
  }else{
    df <- final
    
    #getting the starting coordinate
    myIntervals_start <- df[, c(1,2)] 
    
    #getting the stopping/ending coordinate
    myIntervals_stop <- df[, c(1,3)] 
    
    # Removing "chr" from the chr column 
    myIntervals_start$chr <- myIntervals_start$chr %>% str_remove("chr")
    myIntervals_stop$chr <- myIntervals_stop$chr %>% str_remove("chr")
    
    # Adding additional base to create an arteficial range of coordinates to extract gene names at the locus
    myIntervals_start <- myIntervals_start %>%
      mutate(end = as.numeric(myIntervals_start$start))
    
    myIntervals_stop <- myIntervals_stop %>%
      mutate(end = as.numeric(myIntervals_stop$stop))
    
    # changing colnames
    colnames(myIntervals_start) <- header
    
    colnames(myIntervals_stop) <- header
    
    # The function "makeGRangesFromDataFrame" from the library 
    # GenomicRanges makes an object GRanges called "intervals" from "myIntervals"
    interval_start <- GenomicRanges::makeGRangesFromDataFrame(myIntervals_start)
    intervals_stop = GenomicRanges::makeGRangesFromDataFrame(myIntervals_stop)
    
    # extract the list of all gene_Ids from the txdb object
    genes = genes(txdb)
    
    
    # Use the "annotateIntervals" funtion in order to annotate 
    myAnnotation_start <- as.data.frame(annotateIntervals(interval_start, txdb))
    
    myDf_master_start <- data.frame()
    for (i in 1 :length(myAnnotation_start$gene_id)){
      if(length(c(na.omit(myAnnotation_start$gene_id[i])[[1]])) != 0) {
        myDf_start <- data.frame(chr = myAnnotation_start$seqnames[i], start = myAnnotation_start$start[i],
                                 end = myAnnotation_start$end[i], startAnno = toString(c(na.omit(myAnnotation_start$gene_id[i])[[1]])))
        myDf_master_start <- rbind(myDf_master_start, myDf_start)
      }
      if (length(c(na.omit(myAnnotation_start$gene_id[i])[[1]])) == 0){
        myDF_start_0 <- data.frame(chr = myAnnotation_start$seqnames[i], start = myAnnotation_start$start[i],
                                   end = myAnnotation_start$end[i], startAnno = toString(paste("intergenic")))
        myDf_master_start <- rbind(myDf_master_start, myDF_start_0)
      }
    }
    
    newDF_start <- myDf_master_start[, c(2,4)]
    
    myAnnotation_stop <- as.data.frame(annotateIntervals(intervals_stop, txdb))
    
    myDf_master_stop <- data.frame()
    for (i in 1 :length(myAnnotation_stop$gene_id)){
      if(length(c(na.omit(myAnnotation_stop$gene_id[i])[[1]])) != 0) {
        myDf_stop <- data.frame(chr = myAnnotation_stop$seqnames[i], start = myAnnotation_stop$start[i],
                                end = myAnnotation_stop$end[i], startAnno = toString(c(na.omit(myAnnotation_stop$gene_id[i])[[1]])))
        myDf_master_stop <- rbind(myDf_master_stop, myDf_stop)
      }
      if (length(c(na.omit(myAnnotation_stop$gene_id[i])[[1]])) == 0){
        myDF_stop_0 <- data.frame(chr = myAnnotation_stop$seqnames[i], start = myAnnotation_stop$start[i],
                                  end = myAnnotation_stop$end[i], startAnno = toString(paste("intergenic")))
        myDf_master_stop <- rbind(myDf_master_stop, myDF_stop_0)
      }
    }
    
    newDF_stop <- myDf_master_stop[, c(2,4)]
    newDF_stop %>% dplyr::rename(stop = start,
                                 stopAnno = startAnno) -> newDF_stop
    
    big <- cbind(df, newDF_start)
    big2 <- cbind(big, newDF_stop)
    big3 <- big2 %>% dplyr::select(unique(colnames(.))) %>% relocate(startAnno, stopAnno, .after = stop) %>% as.data.frame()
    
    
    ## write output
    write.csv(big3, file=paste("./results/", id,"_result_final.csv", sep=""), row.names = F)
  }
  
  
  ## final cleanning up without filter out gene of interest
  final1 <- bb4 %>% 
    dplyr::rename(chr=V1,
                  start=V2,
                  stop=V3,
                  readStrand=V4.1,
                  intron_motif=V5,
                  unique_count = V7,
                  multimapped_count=V8,
                  max_overhang = V9,
                  geneStrand=V13,
                  transcriptName=V14,
                  type=V15,
                  totalEx=V16,
                  ExNum=V17,
                  geneName=V18) %>%
    mutate(intron_motif=case_when(
      intron_motif == 0 ~ "non-canonical",
      intron_motif == 1 ~ "GT/AG",
      intron_motif == 2 ~ "CT/AC",
      intron_motif == 3 ~ "GC/AG",
      intron_motif == 4 ~ "CT/GC",
      intron_motif == 5 ~ "AT/AC",
      intron_motif == 6 ~ "GT/AT")) %>% 
    dplyr::select(-V6, -V10, -V11, -V12, -V19, -V20, -V4, 
                  -overlap, -total, -abnormality, -ExNum, -type) %>%
    mutate(sampleName=id) %>%
    as.data.frame()
  
  ## write output
  write.csv(final1, file=paste("./results/", id,"_result_final_allgenes.csv", sep=""), row.names = F)
  
}
