# Parse pattern1 and pattern 2

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use lib qw(~/perl_modules);
use List::MoreUtils qw/ uniq /;

# Prepare input file
# module load tabix
# bgzip -c file.vcf > file.vcf.gz
# tabix -p vcf file.vcf.gz
#  Now use vcftools to make a tab delimited file:
# module load StdEnv/2020 vcftools/0.1.16
# zcat file.vcf.gz | vcf-to-tab > out.tab

# on computecanada:
# module load perl/5.30.2
# module load gcc/9.3.0
# cpan
# install List::MoreUtils

#  This program reads in a tab delimited genotype file generated
#  by vcftools (vcf2tab) and quantifies the number of two pattern types in four taxa

# the four taxa include two diploids (d1, d2) and two tetraploids (t1, t2)

# pattern 1 and pattern 2 required d1 and d2 be homozygous
# pattern 1 and pattern 2 required either t1 or t2 be heterozygous and the other be homozygous
# when t1 is the heterozygous genotype, pattern 1 is defined as d1 being homozygous for one variant in t1 
# and d2 and t2 are homozygous for the other variant
# when t1 is the heterozygous genotype, pattern 2 is defined as d2 being homozygous for one variant in t1 
# and d1 and t2 are homozygous for the other variant
# Rationale: 
# (1) if t1 arose from an ancestor of d1 and d2 before population structure, then pattern 1 and 2 should
# have similar frequencies
# (2) if t1 arose from d2 after d1 and d2 diverged, then pattern 2 should be more abundant than pattern 1
# (3) if t1 and t2 are derived from a common ancestor, pattern 1 should have similar frequencies
# when t1 or t2 are the heterozygous tetraploid, and the same for pattern 2, with differences between pattern 1 and 2
# due to variation between t1 and t2 in the effective population size



# to execute type 
# perl Parse_pattern_1_and_2_tetraploids.pl inputfile.tab 12000000000000000034 pattern12.out 
# where 12000000000000000034 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) d1, (2) d2, (3) t1, (4) t2, or (0) not included

# perl Parse_pattern_1_and_2_tetraploids.pl  combined_Chr7.g.vcf.gz_Chr7_GenotypedSNPs.vcf.gz_filtered.vcf.gz.tab 10200000000000000034 pattern12_10200000000000000034.out


my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $outputfile = $ARGV[2];


unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";
print OUTFILE "d1\td2\tt1\tt2\n";

my @sexes = split("",$ARGV[1]);
my @temp;
my $temp1;
my $temp2;
my $temp3;
my $temp4;
my $number_of_individuals_included=0;
my $pattern1=0;
my $pattern2=0;
my $pairwise_12=0;
my $pairwise_13=0;
my $pairwise_23=0;

my @indiv1=();
my @indiv2=();
my @indiv3=();
my @indiv4=();
my $y;
my @unique_indiv1_nucleotides;
my @unique_indiv2_nucleotides;
my @unique_indiv3_nucleotides;
my @unique_indiv4_nucleotides;

my @display;

my $counter = 0;


for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if(($sexes[$y] == 1)||($sexes[$y] == 2)||($sexes[$y] == 3)||($sexes[$y] == 4)){
		$number_of_individuals_included +=1;
	}	
}	
print "This number should be four: ",$number_of_individuals_included,".\n";
print OUTFILE "d1\td2\tt1\tt2\tpattern1\tpattern2\n";
while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split /[\t\/]/,$line;
	if($temp[0] ne '#CHROM'){
		if($#temp ne (($#sexes+1)*2)+2){
			print "The number of individuals in the input line does not match the number of individuals genotyped ",
			$temp[0],"\t",$temp[1],"\t",$#temp," ",(($#sexes+1)*2)+2,"\n";
		}

		# parse the bases in all genotypes in each sex
		@indiv1=();
		@indiv2=();
		@indiv3=();
		@indiv4=();
		@display=();
		$counter=0;
		for ($y = 3 ; $y <= $#temp; $y=$y+2 ) {
			if(($temp[$y] ne ".")&&($temp[$y+1] ne ".")&&($temp[$y] ne "*")&&($temp[$y+1] ne "*")){ # not allowing asterisks (*), which are deletions
				if($sexes[$counter] == 1){
						push(@indiv1, $temp[$y]);
						push(@indiv1, $temp[$y+1]);
						$temp1=$temp[$y]."/".$temp[$y+1];
				}
				elsif($sexes[$counter] == 2){
					push(@indiv2, $temp[$y]);
					push(@indiv2, $temp[$y+1]);
					$temp2=$temp[$y]."/".$temp[$y+1];
				}
				elsif($sexes[$counter] == 3){
					push(@indiv3, $temp[$y]);
					push(@indiv3, $temp[$y+1]);
					$temp3=$temp[$y]."/".$temp[$y+1];
				}
				elsif($sexes[$counter] == 4){
					push(@indiv4, $temp[$y]);
					push(@indiv4, $temp[$y+1]);
					$temp4=$temp[$y]."/".$temp[$y+1];
				}
			}
			$counter+=1;
		}	
		# OK I should have all the bases loaded for non-missing genotypes for each individual
		
		# find out what and how many unique nucleotides are in each sex
		@unique_indiv1_nucleotides = uniq @indiv1;
		@unique_indiv2_nucleotides = uniq @indiv2;
		@unique_indiv3_nucleotides = uniq @indiv3;
		@unique_indiv4_nucleotides = uniq @indiv4;
		if(($#unique_indiv1_nucleotides == 0)&&($#unique_indiv2_nucleotides == 0)&&
		   ($#unique_indiv3_nucleotides == 1)&&($#unique_indiv4_nucleotides == 0)){
			# this means that d1, d2, and t2 are homozygous and t1 is heterozygous
			if(($unique_indiv1_nucleotides[0] ne $unique_indiv2_nucleotides[0])&&
			($unique_indiv1_nucleotides[0] eq $unique_indiv4_nucleotides[0])&&
			(($unique_indiv2_nucleotides[0] eq $unique_indiv3_nucleotides[0])||
			($unique_indiv2_nucleotides[0] eq $unique_indiv3_nucleotides[1]))
			){
				#print "pattern 2 ",$temp1," ",$temp2," ",$temp3," ",$temp4,"\n";
				$pattern2+=1
			}	
			if(($unique_indiv1_nucleotides[0] ne $unique_indiv2_nucleotides[0])&&
			($unique_indiv2_nucleotides[0] eq $unique_indiv4_nucleotides[0])&&
			(($unique_indiv1_nucleotides[0] eq $unique_indiv3_nucleotides[0])||
			($unique_indiv1_nucleotides[0] eq $unique_indiv3_nucleotides[1]))
			){
				#print "pattern 1 ",$temp1," ",$temp2," ",$temp3," ",$temp4,"\n";
				$pattern1+=1
			}	
		} # end of check that there is at least one genotype in each sex
	} # end of check to see if we are at the first line	
	else{ # print the names of the included samples to the outfile
		for ($y = 0 ; $y <= $#sexes ; $y++ ) {
			if($sexes[$y] == 1){
				$temp1 = $temp[$y+3];
			}	
			elsif($sexes[$y] == 2){
				$temp2 = $temp[$y+3];
			}	
			elsif($sexes[$y] == 3){
				$temp3 = $temp[$y+3];
			}	
			elsif($sexes[$y] == 4){
				$temp4 = $temp[$y+3];
			}	
		}
		#print  $temp1,"\t",$temp2,"\t",$temp3,"\t",$temp4,"\n";
		print OUTFILE $temp1,"\t",$temp2,"\t",$temp3,"\t",$temp4,"\t";
	}
} # end while
print OUTFILE $pattern1,"\t";
print OUTFILE $pattern2,"\n";
close OUTFILE;
```
