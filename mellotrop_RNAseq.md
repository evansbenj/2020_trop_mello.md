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

```
module load StdEnv/2020  gcc/9.3.0
module load kallisto/0.46.1
# index the transcriptome
kallisto index -i /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/mellotrop_RNA/assembly/melltrop_transcriptome_trinityOut
```



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

