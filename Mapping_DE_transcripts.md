# Get fasta seqs from assembly
Working dir:
```
/home/ben/projects/rrg-ben/ben/2020_mellotrop_RNA/mellotrop_RNA/assembly
```

According to this link: https://edwards.sdsu.edu/research/perl-one-liner-to-extract-sequences-by-their-identifer-from-a-fasta-file/ Get fasta entries like this:
```
perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw(42732422)}print if $c' melltrop_transcriptome_trinityOut_3.Trinity.fasta
```
or with a file:
```
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ids.file melltrop_transcriptome_trinityOut_3.Trinity.fasta
```
or
```
awk -v seq="42732422 10004 289064 42567259-,...,42597465-" -v RS='>' '$1 == seq {print RS $0}' melltrop_transcriptome_trinityOut_3.Trinity.fasta
```
