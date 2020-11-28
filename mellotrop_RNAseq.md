# X. mellotrop RNAseq


In this directory we have the raw data:
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/mellotrop_RNA
```
I am trimming it with trimmomatic:
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
LLUMINACLIP:~/projects/rrg-ben/ben/2017_SEAsian_macaques/bin/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TR
AILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  fi
done
```
