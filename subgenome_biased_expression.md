# Assay for subgenome biased expression

* Path
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/mellotrop_RNA/
```

* Index genomes for STAR alignment
```
#!/bin/sh
#SBATCH --job-name=STAR_index
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=256gb
#SBATCH --output=STAR_index.%J.out
#SBATCH --error=STAR_index.%J.err
#SBATCH --account=def-ben

module load StdEnv/2020 star/2.7.9a

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome \
--genomeFastaFiles /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome.fasta \
--sjdbGTFfile /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_Xenbase.gff3 \
--sjdbOverhang 99 \
--limitGenomeGenerateRAM=124544990592

```
* Align genomes with STAR
```
#!/bin/sh
#SBATCH --job-name=STAR_map
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=64gb
#SBATCH --output=STAR_map.%J.out
#SBATCH --error=STAR_map.%J.err
#SBATCH --account=def-ben

module load StdEnv/2020 star/2.7.9a

STAR --genomeDir /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/ \
--runThreadN 6 \
--readFilesIn ${1} ${2} \
--outFileNamePrefix ${3} \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--readFilesCommand zcat

```
* Sort and add readgroups
```
#!/bin/sh
#SBATCH --job-name=sort_and_rg
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=32:00:00
#SBATCH --mem=16gb
#SBATCH --output=sort_and_rg.%J.out
#SBATCH --error=sort_and_rg.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this                                                                                
# sbatch 2022_samtools_sort_and_readgroups.sh bamfile_prefix samplename

module load samtools/1.10
# Sort both alignments
samtools sort ${1}.bam -o ${1}.sorted.bam

# add readgroups
module load picard/2.23.3

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${1}.sorted.bam O=${1}.sorted_rg.bam RGID=4 RGLB=${2} RGPL=ILLUMI
NA RGPU=${2} RGSM=${2}

# index
module load StdEnv/2020 samtools/1.12
samtools index ${1}.sorted_rg.bam

```
* Combine into multibam file with samtools
```
#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:00:00
#SBATCH --mem=8gb
#SBATCH --output=samtools_merge.%J.out
#SBATCH --error=samtools_merge.%J.err
#SBATCH --account=def-ben

module load StdEnv/2020 samtools/1.12

samtools merge -o all_STAR_multibam.bam BJE3795BRAINAligned.sortedByCoord.out.sorted_rg.bam BJE3801BRAINAligned.sortedByCoord
.out.sorted_rg.bam BJE3795HEARTAligned.sortedByCoord.out.sorted_rg.bam BJE3801HEARTAligned.sortedByCoord.out.sorted_rg.bam BJ
E3795LIVERAligned.sortedByCoord.out.sorted_rg.bam BJE3801LIVERAligned.sortedByCoord.out.sorted_rg.bam BJE3796BRAINAligned.sor
tedByCoord.out.sorted_rg.bam BJE3802BRAINAligned.sortedByCoord.out.sorted_rg.bam BJE3796HEARTAligned.sortedByCoord.out.sorted
_rg.bam BJE3802HEARTAligned.sortedByCoord.out.sorted_rg.bam BJE3796LIVERAligned.sortedByCoord.out.sorted_rg.bam BJE3802LIVERA
ligned.sortedByCoord.out.sorted_rg.bam BJE3797BRAINAligned.sortedByCoord.out.sorted_rg.bam BJE3803BRAINAligned.sortedByCoord.
out.sorted_rg.bam BJE3797HEARTAligned.sortedByCoord.out.sorted_rg.bam BJE3803HEARTAligned.sortedByCoord.out.sorted_rg.bam BJE
3797LIVERAligned.sortedByCoord.out.sorted_rg.bam BJE3803LIVERAligned.sortedByCoord.out.sorted_rg.bam BJE3800LIVERAligned.sort
edByCoord.out.sorted_rg.bam
```
