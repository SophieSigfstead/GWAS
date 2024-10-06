# GWAS
** This is a WIP and will be completed on October 6*

This document provides some details on the repository. This repository contains many trials of our GWAS SNP filtering method.All files associated with a specific gwas are labelled beginning with “gwas_#”. To reference which gwas study this corresponds to please look at references.txt.

The main protocol for each GWAS analysis is to:
1. Identify an Enformer track list relevant to the trait being studies
2. Merge 1000 genomes SAD scores for each track in the track list with summary statistics for the given GWAS using a gwas_#_matching.py. These files match rows of summary statistics by either snp id (rs id) or chromosome and position pairs, depending on the information provided in the summary statistics.This is done using a matching routine. For gwas 1, this is called "gwas_1_matching.py". For gwas 3, this is called "gwas_3_scz_matching.py". For gwas 4, this is called "gwas_4_matching.py".
3. RReplicate the given GWAS using a pre-filtered list of SNPs. These files are labelled gwas_#_leading_snps_procedure.py. These files first filter the SNPs based upon global thresholds for Enformer SAD track scores, and then adjust the p-value based upon the number of SNPs filtered from the given study. This is done with gwas_1_leading_snps_procedure.py for gwas 1, and gwas_3_scz_leading_snps_procedure.py for gwas 3, and gwas_4_alz_leading_snps_procedure.py for gwas_4.





