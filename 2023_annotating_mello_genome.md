# Annotating a new mello assembly using trop data
Path on graham:
```
/home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome
```
First extract exons from trop gff file:
```
grep 'exon' XENTR_10.0_Xenbase.gff3 > XENTR_10.0_Xenbase_exonsonly.gff3
```
Now use this to extract fasta seqs for each exon:
```
module load bedtools
bedtools getfasta -fi XENTR_10.0_genome.fasta -bed XENTR_10.0_Xenbase_exonsonly.gff3 -fo XENTR_10.0_genome_exonsonly.fasta
```
This yeilds 648717 exons from trop.

Now use this to identify top hits for each trop exon:
```
module load StdEnv/2020 minimap2/2.24
minimap2 -t8 reference.fasta query.fasta >alignments.sam 
```

THis is the mello genome assembly I am working with (for now at least):
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/assembly_Germanydata_10kb_matepair_plus_pacbio_genome.fa
```
