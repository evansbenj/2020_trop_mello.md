# Masurca assembly

path:
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/assembly_Germanydata_10kb_matepair_plus_pacbio/CA.mr.97.17.15.0.02
``

finished genome assembly using Germanydata_10kb_matepair_plus_pacbio:
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/assembly_Germanydata_10kb_matepair_plus_pacbio_genome.fa
```

This genome assembly only has the Germany paired end data plus the pacbio seqs. It did not use the supernova assembly; those efforts failed on graham and beluga and never got started on cedar


```
#!/bin/sh 
#SBATCH --job-name=masurca
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=32
#SBATCH --time=168:00:00 
#SBATCH --mem=500gb
#SBATCH --output=masurca.%J.out
#SBATCH --error=masurca.%J.err
#SBATCH --account=def-ben

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load StdEnv/2020  gcc/9.3.0 masurca/4.1.0
./assemble.sh
```
