# Annotating a new mello assembly using trop data
Path on graham:
```
/home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome
```
First extract exons from trop gff file:
```
grep 'exon' XENTR_10.0_Xenbase_longest.gff3 > XENTR_10.0_Xenbase_longest_exonsonly.gff3
```
Now make a new bed file that also has the name of each exon in it:
```
cut -f1,4,5,9 XENTR_10.0_Xenbase_longest_exonsonly.gff3 > XENTR_10.0_Xenbase_longest_exonsonly_names.gff
```
Now use this to extract fasta seqs for each exon:
```
module load bedtools
bedtools getfasta -name -fi XENTR_10.0_genome.fasta -bed XENTR_10.0_Xenbase_longest_exonsonly_names.gff -fo XENTR_10.0_genome_exonsonly.fasta
```
This yeilds 236189 exons from trop.

Now use this to identify top hits for each trop exon:
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

module load StdEnv/2020 minimap2/2.24

#/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/bin/minimap2/
minimap2 -x asm10 -P /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/assembly_Germanydata_10kb_matepair_plus_pacb
io_genome.fa ${1}.fasta > ${1}.paf
```

This is the mello genome assembly I am working with (for now at least):
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/assembly_Germanydata_10kb_matepair_plus_pacbio_genome.fa
```

# Trying with BLAST

makeblastdb:
```
```
blast long exons against mello assembly:
```
blastn -task dc-megablast -query ../../2020_XT_v10_refgenome/XENTR_10.0_genome_exonsonly.fasta -db assembly_Germanydata_10kb_matepair_plus_pacbio_genome.fa_blastable -outfmt 6 -out XENTR_10.0_genome_exonsonly_to_mello_blast.out
```
Save top two alignments based on bit score (still working on this):
```
cat XENTR_10.0_genome_exonsonly_to_mello_blast.out | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > top_two.out
```

