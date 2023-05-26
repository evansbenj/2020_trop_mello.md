# DBG2OLC

directory:
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/bin/DBG2OLC
```

This assembler works on a single thread and requires a fragmented assembly that is generated from short reads. I used the Germany supernovo assembly for this and added the pacbio for the hybrid assembly.  It worked surprisingly quickly (~2 days).

```
#!/bin/sh 
#SBATCH --job-name=DBG2OLC
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00 
#SBATCH --mem=512gb
#SBATCH --output=DBG2OLC.%J.out
#SBATCH --error=DBG2OLC.%J.err
#SBATCH --account=def-ben

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/bin/DBG2OLC/DBG2OLC/compiled/DBG2OLC k 17 AdaptiveTh 0.005 KmerCovTh 
2 MinOverlap 20 RemoveChimera 1 Contigs /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/Germany_genome/Super_NovaXeno
_mega_gt200.fasta f /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/2023_Mello_PacBio/2017_mellotropicalis_PacBio/all
_PacBio_combined.fasta
```
and part 2:
```
#!/bin/sh 
#SBATCH --job-name=DBG2OLC
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00 
#SBATCH --mem=512gb
#SBATCH --output=DBG2OLC.%J.out
#SBATCH --error=DBG2OLC.%J.err
#SBATCH --account=def-ben
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load nixpkgs/16.09 intel/2018.3 blasr/5.3.0

#this is to concatenate the contigs and the raw reads for consensus

# cat /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/Germany_genome/Super_NovaXeno_mega_gt200.fasta /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/2023_Mello_PacBio/2017_mellotropicalis_PacBio/all_PacBio_combined.fasta > ctg_pb.fasta

# we need to open a lot of files to distribute the above file into lots of smaller files
# ulimit -n unlimited # this does not work 

#run the consensus scripts

sh /home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/bin/DBG2OLC/Sparc/utility/split_and_run_sparc.sh backbone_raw.fasta DBG2OLC_Consensus_info.txt ctg_pb.fasta ./consensus_dir 2 >cns_log.txt

```
