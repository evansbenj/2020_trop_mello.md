# Mapping contigs to XT using minimap2

```
#!/bin/sh
#SBATCH --job-name=minmap
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --output=minmap.%J.out
#SBATCH --error=minmap.%J.err
#SBATCH --account=def-ben

# minimap2 -x asm10 -a --secondary=no -t8 reference.fasta query.fasta >alignments.sam
module load StdEnv/2020 minimap2/2.24

#/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/bin/minimap2/
minimap2 -x asm10 -a --secondary=no -t8 /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome.
fasta /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/assembly_Germanydata_10kb_matepair_plus_pacbio_genome.f
a >alignments.sam 
```

# convert to bam, sort and index
```
#!/bin/sh
#SBATCH --job-name=samtools
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=12gb
#SBATCH --output=samtools.%J.out
#SBATCH --error=samtools.%J.err
#SBATCH --account=def-ben
module load samtools
samtools view -S -b alignments.sam > mello_masurca_to_XTv10.bam
samtools sort mello_masurca_to_XTv10.bam -o mello_masurca_to_XTv10_sorted.bam
samtools index mello_masurca_to_XTv10_sorted.bam
```

# Extract Chr7 and index
```
samtools view -b mello_masurca_to_XTv10_sorted.bam Chr7 > Chr7_mello_masurca_to_XTv10_sorted.bam
samtools index Chr7_mello_masurca_to_XTv10_sorted.bam
```
