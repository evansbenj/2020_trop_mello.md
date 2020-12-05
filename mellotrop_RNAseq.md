# X. mellotrop RNAseq


In this directory we have the raw data:
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/mellotrop_RNA
```
after copying the data, I had to change the permissions of each sample directory to avoid getting a disk quota error like this:
```
chmod g+s Sample_BJE3796BRAIN
```

I am trimming it with trimmomatic (after copying this file to the same directory as the sbatch file `~/projects/rrg-ben/ben/2017_SEAsian_macaques/bin/Trimmomatic-0.36/adapters/TruSeq2-PE.fa`:
```
#!/bin/sh
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=8gb
#SBATCH --output=trimmomatic.%J.out
#SBATCH --error=trimmomatic.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2020_trimmomatic.sh ../raw_data/plate1

module load StdEnv/2020
module load trimmomatic/0.39

#_R1_001.fastq.gz
for file in $1/*_R1*.fastq.gz ; do         # Use ./* ... NEVER bare *
  if [ -e "$file" ] ; then   # Check whether file exists.
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE ${file::-14}1${file:(-13)} ${file::-14}2${file:(-13)} ${file
::-14}1.trimmed.fq.gz ${file::-14}1.trimmed_single.fq.gz ${file::-14}2.trimmed.fq.gz ${file::-14}2.trimmed_single.fq.gz I
LLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:3 TR
AILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  fi
done
```

After trimming all the files, I concatenated the forward reads. And then I concatenated the reverse reads:
```
cat ../S*/*_R1.trimmed_00*.fastq.gz > all_R1.fq.gz
cat ../S*/*_R2.trimmed_00*.fastq.gz > all_R2.fq.gz
```
in this directory:
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/mellotrop_RNA/assembly
```
# Assemble with Trinity


```
#!/bin/sh
#SBATCH --job-name=trinity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00
#SBATCH --mem=256gb
#SBATCH --output=trinity.%J.out
#SBATCH --error=trinity.%J.err
#SBATCH --account=def-ben


module load nixpkgs/16.09  gcc/7.3.0 nixpkgs/16.09
module load openmpi/3.1.2
module load samtools/1.9
module load salmon/0.11.3
module load bowtie2/2.3.4.3
module load jellyfish/2.2.6
module load trinity/2.8.4

Trinity --seqType fq --left /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/mellotrop_RNA/assembly/all_R1.fq.gz --right /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/mellotrop_RNA/assembly/all_R2.fq.gz --CPU 20 --full_cleanup --max_memory 200G --min_kmer_cov 2 --include_supertranscripts --output /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/mellotrop_RNA/assembly/melltrop_transcriptome_trinityOut
```

# Count with Kallisto

Descriptions are here: https://pachterlab.github.io/kallisto/starting.html

```
module load StdEnv/2020  gcc/9.3.0
module load kallisto/0.46.1
# index the transcriptome can be done directly with kallisto (the -i flag tells kallisto what to name the index)
# and the fasta follows this with no flag)

kallisto index -i tropicalis_transcriptome_trinityOut.Trinity.fasta.kallisto_idx ./tropicalis_transcriptome_trinityOut.Trinity.fasta

# or via a perl script that comes with trinity (on info):
perl /usr/local/trinity/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts tropicalis_transcriptome_trinityOut.Trinity.fasta --seqType fa --samples_file samplefile.tsv --est_method kallisto --output_dir ./kallisto_denovo/ --trinity_mode --prep_reference

```
# Now for each file count abundances of each transcript
```
for file in ../../data/trimmed_RNAseq_data/X*_R1_paired.fastq.gz ; do      
  if [ -e "$file" ] ; then   # Check whether file exists.
      r1="${file::${#file}-19}_R1_paired.fastq.gz"
      r2="${file::${#file}-19}_R2_paired.fastq.gz"
      kallisto quant -i ./tropicalis_transcriptome_trinityOut.Trinity.fasta.kallisto_idx -o 
./counts/${file:31:4} ${r1} ${r2}
  fi
done
```

# I think this does the same thing but runs via trinity
```
for file in ../../data/trimmed_RNAseq_data/X*_R1_paired.fastq.gz ; do      
  if [ -e "$file" ] ; then   # Check whether file exists.
      r1="${file::${#file}-19}_R1_paired.fastq.gz"
      r2="${file::${#file}-19}_R2_paired.fastq.gz"
      perl /usr/local/trinity/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance
.pl --transcripts tropicalis_transcriptome_trinityOut.Trinity.fasta --seqType fq --left ${r1} --right ${r2} --est_method kallist
o --aln_method bowtie --samples_file troptad_samples.txt --trinity_mode --output_dir counts_
trinity
  fi
done
```

In the above command, the sample file is this:
```

#      --samples_file <string>    tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
#                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
#                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
#                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
```

```
female	female_rep1	../../data/trimmed_RNAseq_data/XT2_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT2_R2_paired.fastq.gz
female	female_rep2	../../data/trimmed_RNAseq_data/XT3_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT3_R2_paired.fastq.gz
female	female_rep3	../../data/trimmed_RNAseq_data/XT6_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT6_R2_paired.fastq.gz
female	female_rep4	../../data/trimmed_RNAseq_data/XT9_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT9_R2_paired.fastq.gz
female	female_rep5	../../data/trimmed_RNAseq_data/XT10_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT10_R2_paired.fastq.gz
female	female_rep6	../../data/trimmed_RNAseq_data/XT11_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT11_R2_paired.fastq.gz
female	female_rep7	../../data/trimmed_RNAseq_data/XT16_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT16_R2_paired.fastq.gz
female	female_rep8	../../data/trimmed_RNAseq_data/XT17_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT17_R2_paired.fastq.gz
female	female_rep9	../../data/trimmed_RNAseq_data/XT20_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT20_R2_paired.fastq.gz
male	male_rep1	../../data/trimmed_RNAseq_data/XT3_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT3_R2_paired.fastq.gz
male	male_rep2	../../data/trimmed_RNAseq_data/XT9_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT9_R2_paired.fastq.gz
male	male_rep3	../../data/trimmed_RNAseq_data/XT11_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT11_R2_paired.fastq.gz
male	male_rep4	../../data/trimmed_RNAseq_data/XT20_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT20_R2_paired.fastq.gz
male	male_rep5	../../data/trimmed_RNAseq_data/XT7_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT7_R2_paired.fastq.gz
male	male_rep6	../../data/trimmed_RNAseq_data/XT8_R1_paired.fastq.gz	../../data/t
rimmed_RNAseq_data/XT8_R2_paired.fastq.gz
```

# DE analysis with edgeR


#Transcript abundance quantification (from Xue)
```
for i in XXX/*_R1_paired.fastq.gz; 
do name=$(grep -o "XT[0-9]*" < (echo $i));r1=/home/xue/tropicalis_gonad_transcriptome_Dec2018/data/trim/$name\_R1_paired.fastq.gz;r2=/home/xue/tropicalis_gonad_transcriptome_Dec2018/data/trim/$name\_R2_paired.fastq.gz; 
kallisto quant -i /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/transcript_expression_raw_count/kallisto_indexing_transcriptome/tropicalis_transcriptome_trinityOut.Trinity.fasta.kallisto_idx  -o /home/xue/tropicalis_gonad_transcriptome_Dec2018/analysis/transcript_expression_raw_count//$name <(gunzip -c $r1) <(gunzip -c $r2);done
```
#To compute the matrix (time cost: ~10min)-didnt do yet
```
time perl /home/xue/software/trinityrnaseq-Trinity-v2.5.1/util/abundance_estimates_to_matrix.pl --est_method kallisto --out_prefix tropicalis_gonad --gene_trans_map /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut/tropicalis_gonad_supertranscriptome_dec2018/tropicalis_transcriptome_trinityOut.Trinity.fasta.gene_trans_map --name_sample_by_basedir XT1/abundance.tsv XT10/abundance.tsv XT11/abundance.tsv XT13/abundance.tsv XT16/abundance.tsv XT17/abundance.tsv XT19/abundance.tsv XT2/abundance.tsv XT20/abundance.tsv XT3/abundance.tsv XT6/abundance.tsv XT7/abundance.tsv XT8/abundance.tsv XT9/abundance.tsv
```


Now quantify for each sample:
```
kallisto quant -i /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/mellotrop_RNA/assembly/melltrop_transcriptome_trinityOut.idx -o Sample_name XX.R1.fq.gz XX.R2.fq.gz

# -o Sample_name: Sample_name will be name of the output folder produced from this analysis
# XX.R1.fq.gz XX.R2.fq.gz : input FASTQ files
```


https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

Below is the R script from Xue that was used for edgeR analysis of the trop tads.  I need to modify this so that we examine the interaction between sex biased expression and tissue type (brain, liver, heart).
```
setwd("/Users/Shared/Previously Relocated Items/Security/benstuff/  teaching/BIO720_2019/RNAseq_example")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library(tidyverse)
library(edgeR)

#read in the data

#build a string that defines the path and file name of the count data
input_file_folder <- file.path(".")
input_file <- paste(input_file_folder,"tropicalis_gonad.isoform.counts.matrix", sep="/")

# read the count data into a dataframe called "express_data"
express_data <- read.table(input_file, header=T, row.names = 1, com='')

# get rid of any rows that have incomplete data
express_data <- express_data[complete.cases(express_data), ]

#Define subsets of individuals for analysis
# all individuals
stage50 <- c("XT2","XT3","XT6","XT9", "XT10","XT11","XT16","XT17", "XT20","XT1","XT7","XT8","XT13","XT19" )
# WW females, WY males
group1 <- c("XT3", "XT9", "XT11","XT20", "XT7", "XT8")

# WZ females, ZY males
group2 <- c("XT2", "XT6", "XT10", "XT16", "XT17", "XT1", "XT13", "XT19")

# WW females, ZY males
group3 <- c("XT3", "XT9", "XT11","XT20","XT1", "XT13", "XT19")

# WZ females, WY males
group4 <- c("XT2", "XT6","XT10", "XT16", "XT17","XT7", "XT8")

# Define sex for each analysis
stage50_sex <- factor(c(rep("female", 9), rep("male", 5)))
group1_sex <- factor(c(rep("female", 4), rep("male", 2)))
group2_sex <- factor(c(rep("female", 5), rep("male", 3)))
group3_sex <- factor(c(rep("female", 4), rep("male", 3)))
group4_sex <- factor(c(rep("female", 5), rep("male", 2)))

length(stage)

#select stage and conditions
stage = group4
conditions = group4_sex

# make a matrix for analysis that has only the subset of the data
# that we want to analyze
rnaseqMatrix <- express_data[,match(stage, names(express_data))]
dim(rnaseqMatrix)

# get rid of rows where there is less than an average of one read per individual
# in the analysis
rnaseqMatrix = rnaseqMatrix[rowSums(rnaseqMatrix)> length(stage),]
dim(rnaseqMatrix)

# do DE with edgeR
# make a DGEList object using the count data and the group information 
# (here which sex the sample is)
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)

# calculate normalization factors to scale each
# library size by
exp_study = calcNormFactors(exp_study)

# estimate negative binomial dispersion across
# all genes
exp_study = estimateCommonDisp(exp_study)

# estimate dispersion across
# individual genes as well
exp_study = estimateTagwiseDisp(exp_study)

# perform an exact test for differences between two groups
# of negative binomial counts
et = exactTest(exp_study, pair=c("female", "male"))

# this converts the object "et" into a list
# which is easier to parse. Can also be used to 
# obtain a subset of the loci
tTags = topTags(et, n=NULL, sort.by = "none")
is.list(tTags)

# get the logFC, logCPM, PValue, and FDR from this list
result_table = tTags$table
is.list(result_table) # results table is a list

# convert it to a dataframe and add columns with the treatments defined
result_table = data.frame(sampleA="female", sampleB="male", result_table)

#exp_study$samples$group

# write the output
output_folder <- file.path(".")
output_file <- paste(output_folder,"de_result_trop_uncollapsed_group4_edgeR.tsv", sep="/")
write.table(result_table, output_file, na = " ", sep = "\t", quote = FALSE, row.names = TRUE,col.names = TRUE)

```
