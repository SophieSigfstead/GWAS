# GWAS

This document provides some details on the repository.

The main protocol for each GWAS analysis is to:
1. Identify an Enformer track list relevant to the trait being studies
2. Merge 1000 genomes SAD scores for each track in the track list with summary statistics from the given GWAS. This is done using a matching routine. For gwas 1, this is called "gwas_1_matching.py". For gwas 3, this is called "gwas_3_scz_matching.py".
3. Replicate the GWAS. This is done with gwas_1_leading_snps_procedure.py for gwas 1, and gwas_3_scz_leading_snps_procedure.py for gwas 3. 

 
