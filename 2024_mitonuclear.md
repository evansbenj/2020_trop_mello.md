
# Mitonuclear interactions and subgenome evolution

* get the coding sequence from xenbase
* make sure there is a hard return after the sequence
* make sure it begins with a start (ATG) and a stop codon (TAA, TAG, TGA)

Here is a script to convert fasta files to phylip format (modified from https://github.com/josephhughes/Sequence-manipulation/blob/master/Fasta2Phylip.pl)

```perl
#!/usr/bin/perl -w
# obtained from Yu-Wei's Bioinformatics playground 
# http://yuweibioinfo.blogspot.com/2009/01/fasta-to-phylip-converter.html

use strict;

MFAtoPHYLIP($ARGV[0]);

sub MFAtoPHYLIP
{
	my $inline;
	my $outfile = "$_[0]\.phy";
	my $count = 0;
	my $len;
	my $substate = 0;
	my @subheader;
	my @subcontent;
	my $m;
	my $n;

	open (FILE, "<$_[0]");
	while (defined($inline = <FILE>))
	{
		chomp($inline);
		if ($inline =~ /^>([A-Za-z0-9.\-_:]+)/)
		{
			$subheader[$count] = $1;
			$subcontent[$count] = "";
			$count++;
		}
		else
		{
			$subcontent[$count - 1] = $subcontent[$count - 1] . " $inline";
		}
	}
	close (FILE);

	# Calculate the content length
	$n = length($subcontent[0]);
	$len = $n;
	for ($m = 0; $m < $n; $m++)
	{
		if (substr($subcontent[0], $m, 1) eq " ")
		{
			$len--;
		}
	}

	open (FILE, ">$outfile");
	print FILE "   $count    $len\n";
	for ($m = 0; $m < $count; $m++)
	{
		$len = 10 - length($subheader[$m]);
		my $temp = substr($subheader[$m], 0, 20);
		print FILE "$temp";
		for ($n = 0; $n < $len; $n++)
		{
			print FILE " ";
		}
		print FILE " $subcontent[$m]\n";
	}
	close (FILE);
}
```

# blast XT gene to germany genome
```
module load  StdEnv/2020  gcc/9.3.0 blast+/2.14.0
blastn -query XT_pars2.fa -db /project/6019307/ben/2020_mellotrop_RNA/Germany_genome/Super_NovaXeno_mega_gt200.fasta_blastable -outfmt 6 -out XT_pars2_to_Germany.out
```
# get coordinates of the matches
```
cut -f2,9,10 XT_pars2_to_Germany.out > hitz.bed
```
# modify the hitz.bed to make sure start < end

# get fasta
```
module load bedtools
bedtools getfasta -fi /project/6019307/ben/2020_mellotrop_RNA/Germany_genome/Super_NovaXeno_mega_gt200.fasta -bed hitz.bed -fo mel_pars2.fa
```

# align 
```
module load StdEnv/2020 mafft/7.471
cat XT_pars2.fa mel_pars2.fa > all_pars2.fa
```
# convert to phylip format
```
./fasta2phylip.pl all_pars2_aligned.fa
```

# get rid of periods in the names
```
sed -i -e 's/\./_/g' all_pars2_aligned.fa.phy
```

# check distance
```
module load StdEnv/2020  intel/2020.1.217 phylip/3.698
dnadist
```
