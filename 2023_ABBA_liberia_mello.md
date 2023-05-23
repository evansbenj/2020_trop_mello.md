# Testing Liberia species status using site patterns

For a genomic position where there is a divergent site between each mellotrop homeolog, there are three possible (interesting) site patterns that correspond to different scenarios.
* If trop + liberia are the same species, then the phylogeny is (((T,L),M1),M2) and the site pattern (in the order T,L,M1,M2) ABBA should be the same as BABA, where B is derived and A is ancestral
* If trop is more closely related to M1, then the phylogeny is (((T,M1),L),M2) and the site pattern (in the order T,L,M1,M2) BABA should more frequent than ABBA, where B is derived and A is ancestral
* If liberia is more closely related to M1, then the phylogeny is (((L,M1),T),M2) and the site pattern (in the order T,L,M1,M2) ABBA should more frequent than BABA, where B is derived and A is ancestral

This can be tested by 
- searching for positions that are heterozygous in the mellotrop genotype and then 
- searching within these sites for positions that are also different between trop and liberia
- counting patterns and comparing these counts
