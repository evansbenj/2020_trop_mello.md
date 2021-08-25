# Goals and questions for mellotrop project
* What is the level of divergence between the subgenomes?
* Is there a different population of TEs in each subgenome?
* Is there a different rate of TE mobility in each subgenome?
* Is there a different rate of gene silencing in each subgenome?
* Is there a different level of expression in each subgenome (this is difficult to assess directly)?
* Do pseudogenized genes in XL also tend to be pseudogenized in X. mellotrop?

## Step 1: For each contig, estimate distance to trop genome
#### Expectation: there should be a bimodal distribution of distances that correspond to each subgenome.  Depending on the level of divergence, these distributions probably will overlap
#### Methods: possibly Blast, gmap, bwa???

## Step 2: For each genomic region in trop, assess how many contigs map to that region; use distances to infer which subgenome they are derived from
#### Challenges include chimerical contigs (that are made of sequence from each subgenome), also heterozygosity could cause issues
#### work from results from Step 1

