# ComputationalBiology
Original code scripted for Computational Biology course assigments. 
Python .py code written and de-bugged using PyCharm, and subsequently written up as a Jupyter Notebook.
R .rmd code written and tested in Markdown format using R Notebook in RStudio.

Assigment1:
Python functions.
1) For a given DNA sequence, calculates the frequency of each DNA base.
2) For a given DNA sequence, calculates the frequency of a given base pair.
3) For a given DNA sequence, calculates the frequency of ALL possible base pairs.

Assignment2:
Python functions.
1) For a given .fastq next generation squencing file, for each posiiton of the sequence, calculates the fraction of reads with quality scores greater than or equal to 30 at that posision.
2) For a given .fastq next generation squencing file, for all possible k's, calculates the number of reads with exactly k positions with quality scores greater than or equal to 30.

Assignment3:
Python functions.
For a given sequence string, finds all open reading frames (ORF) longer than a given threshold value.
ORF length defined as the number of codons (3 bases), including the start ("ATG") and stop ("TAG", "TAA", "TGA") codons.

Assignment4:
Python functions.
For a given multiple sequence alignment file, calculates the number of positions where all 4 alleles (A, C, G, T) appear.

Assignment5:
R functions.
1) For a given file from the '1000 Genomes Project', calculates the number of single nucleotide polymorphisms (SNPs) for which the alternative allele is present for a given sample of individuals.
Calculated this number for the sample of European individuals.
Calculated this number for the sample of African individuals.
2) For a given file from the '1000 Genomes Project', computes the allele frequency spectrum for a given sample of individuals.
For a given sample, the allele frequency spectrum details the number of SNPs for which there are k haplotypes that have the alternative allele.
Computed and plotted the allele frequency spectrum for the sample of European individuals and for the sample of African individuals.

Assignment6:
R functions:
1) For a given file from the '1000 Genomes Project', calculates the average pairwise diversity between two samples of individuals.
The pairwise diversity number details how many SNPs differ between the haplotypes of a random individual from sample1 and a random indivividual from sample2.
The mean pairwise diversity is the average of the mean of the differences for n randomly chosen pairs.
2) For a given file from the '1000 Genomes Project', computes the 95% confidence interval for the average pairwise diversity between two samples of individuals.
Calculated this within the sample of European individuals.
Calculated this within the sample of African individuals.
Calculated this between the sample of European and the sample of African individuals.

Referenced analyzed files also included:
XI1_ATCACG_L001_R1_001.fastq,
XI1_ATCACG_L001_R2_001.fastq,
RETT-1_S1_L001_R1_001.fastq,
RETT-1_S1_L001_R2_001.fastq,
ACE2.fasta,
part-ace2-multi.txt,
abbgen1k.csv,
shorttesteg.csv

Source file for referenced function:
lecture2functions.py