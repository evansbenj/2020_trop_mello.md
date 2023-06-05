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
minimap2 -x asm10 -a --secondary=no -t8 reference.fasta query.fasta >alignments.sam 
```
