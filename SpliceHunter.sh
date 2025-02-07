#!/bin/bash

################################################################################################################
## Filename: SpliceChaser.sh
## Created: October 10, 2024
## Author(s): Muneeza Maqsood, Chung Kok
##
## Purpose: 
##      - to execute SpliceChaser Pipeline with conversion format implemented
##          
## Instructions:  
##      - put all SJ.out.tab files in the input directory
##      - output atypical splicing events
##
## Time estimation:
##      - 
##      - 
##         
################################################################################################################


############ parameter to set #################################
## specify output parent and tmp directory
if [ ! -d "./input/tmp" ] 
then
    mkdir -pv ./input/tmp
    echo "created tmp folder in input directory"
fi

if [ ! -d "./results" ] 
then
    mkdir -pv ./results
    echo "created results directory"
fi

if [ ! -d "./final_results" ] 
then
    mkdir -pv ./final_results
    echo "created final results directory"
fi

## specify bams and tmp_bams directory containing bams and megadeoth's output
if [ ! -d "./bams/tmp_bams" ] 
then
    mkdir -pv ./bams/tmp_bams
    echo "created tmp_bams folder in bams directory"
fi
## specify output directory to store shannon_diversity output   
if [ ! -d "./bams/tmp_bams/output" ] 
then
    mkdir -pv ./bams/tmp_bams/output
    echo "created output folder in tmp_bams directory"
fi
if [ ! -d "./bams/tmp_bams/output/NovelSplices" ] 
then
    mkdir -pv ./bams/tmp_bams/output/NovelSplices
    echo "created NovelSplices folder in output directory"
fi

## pre-determined files and folders
inputDir=./input
dataDir=./input/tmp
prep=./script/pulling_true_pos_from_megadepth.R
db=../../database/RefSeq_ref_with_Intron_v1.bed
cf=./script/convertFileFormat.R
stage1=./script/junction_extraction_RefSeq.R
stage2=./script/PopulationSummarySJ.R
bamsDir=./bams
tmpDir=./bams/tmp_bams
md=../../database/tools/megadepth

## converting crams to bams
cd $bamsDir

echo -e "\n convert crams to bams... \n"

if [ $(ls | grep -c ".cram$") == 0 ]; then echo -e "\n No crams found \n"; else for i in $(ls *.cram); do id=${i%%.cram}; echo $id; samtools view -b -T ../database/WholeGenomeFasta/genome.fasta -o ./tmp_bams/$id.bam $i; done; fi
cd ../

cd $tmpDir

echo -e "\n Running Megadepth... \n"

if [ $(ls | grep -c "jxs.tsv$") == 0 ]; then for i in $(ls *.b37.bam); do id=${i%%.cram}; echo $id; $md $i --all-junctions; done; else echo -e "\n Megadepth output already exists \n"; fi

cd ../..

## Pulling variants (small deletions of upto 50bp) from MegaDepth output (if not already present in STAR's)
echo -e "\n Extraction of missing deletion variants ... begin \n "
Rscript --vanilla $prep

## convert file format first
echo -e " Conversion format process ... begin \n "
Rscript --vanilla $cf

## start running the program
cd $dataDir

echo -e "\n SpliceChaser program initiating ... \n "
for i in $(ls *.bed)
do
    # extract filename
    id=${i%%.bed}
	
	echo "processing ... $id"

	## perform bedtools
	bedtools intersect -a $i -b $db -wo > X__${id}_mapped.bed
done

## back to parent directory
cd ../..

## stage 1
echo -e "\n Stage 1 initiating ... \n "
Rscript --vanilla $stage1


## removing tmp folder
echo -e " Files cleanning ... \n "
rm -rf $dataDir

## Calculating Shannon_diversity Index

echo -e "\n Shannon Diversity Calculations Begin... \n"
python3.10 ./script/Shannon_Score_Calc.py

## stage 2
echo -e " Stage 2 initiating - Population Summary ... \n "
Rscript --vanilla $stage2

echo -e "\n Analysis ... completed"

