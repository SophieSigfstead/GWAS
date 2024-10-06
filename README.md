# GWAS Repository Information

This document provides some details on the repository. This repository contains many trials of our GWAS SNP filtering method. 

The main protocol for each GWAS analysis is to:
1. Identify an Enformer track list relevant to the trait being studied. Tracks are labelled according to tissue type in target_dnase_ataq_tracks_labelled.csv.

2. Merge 1000 genomes SAD scores for each track in the track list with summary statistics for the given GWAS using a gwas_#_matching.py. These files match rows of summary statistics by either snp id (rs id) or chromosome and position pairs, depending on the information provided in the summary statistics. For gwas 1, this is called "gwas_1_matching.py". For gwas 3, this is called "gwas_3_scz_matching.py". For gwas 4, this is called "gwas_4_matching.py".

3. Replicate the given GWAS using a pre-filtered list of SNPs. These files are labelled gwas_#_leading_snps_procedure.py. These files first filter the SNPs based upon global thresholds for Enformer SAD track scores, and then adjust the p-value based upon the number of SNPs filtered from the given study. This is done with gwas_1_leading_snps_procedure.py for gwas 1, and gwas_3_scz_leading_snps_procedure.py for gwas 3, and gwas_4_alz_leading_snps_procedure.py for gwas 4.

Note that all files associated with a specific gwas are labelled beginning with “gwas_#”. To reference which gwas study this corresponds to please look at references.txt. Note that the analysis for gwas 2 was not completed due to limitations in study replicability. The following analyses have been conducted, as of Oct 6, 2024:

### GWAS 1:


Results of each threshold t are stored in GWAS/filtered_snps_gwas_1/filtered_snps_gwas_1_threshold=<t>.csv.

Analysis of results using the following data:
- eQTL (Data available in GWAS/EQTLs)
- FANTOM5_annotations (Data available in GWAS/FANTOM5_annotations)

### GWAS 2:

Results of each threshold t are stored in GWAS/filtered_snps_gwas_2/filtered_snps_gwas_2_threshold=<t>.csv.

### GWAS 3:

### GWAS 4:

Other Directories:
### GWAS/enformer_and_tangermeme 

### Finemapping: GWAS/finemapping and GWAS/fine_mapping_results

### Genome Wide Complex Trait Analysis (GCTA): GWAS/gcta_output_directory and GWAS/gcta
- GCTA software website: https://yanglab.westlake.edu.cn/software/gcta/#Overview
- GCTA install is in GWAS/gcta
- GWAS/gcta_output_directory housed intermediate gcta files while gwas_2 was running, but is currently empty. GWAS 2 was not able to be replicated well, so the analysis was not completed. 

### gencode_annotations

### genome_assembly
- Contains genome annotations. Version v46lift37 (GWAS/genome_assembly/gencode.v46lift37.annotation.gff3) is most applicable to hg37, which is applicable to gwas 1-4.
- Contains make_coding_regions_list.py. This was used to create the exon_regions that are used to identify coding snps in the GWAS filtering procedure. This routine creates the regions for coding snps, which are contained in GWAS/genome_assembly/exon_regions_v2.csv


The main protocol for each GWAS analysis is to:
1. Identify an Enformer track list relevant to the trait being studied. Tracks are labelled according to tissue type in target_dnase_ataq_tracks_labelled.csv.

2. Merge 1000 genomes SAD scores for each track in the track list with summary statistics for the given GWAS using a gwas_#_matching.py. These files match rows of summary statistics by either snp id (rs id) or chromosome and position pairs, depending on the information provided in the summary statistics. This is done using a matching routine. For gwas 1, this is called "gwas_1_matching.py". For gwas 3, this is called "gwas_3_scz_matching.py". For gwas 4, this is called "gwas_4_matching.py".

3. Replicate the given GWAS using a pre-filtered list of SNPs. These files are labelled gwas_#_leading_snps_procedure.py. These files first filter the SNPs based upon global thresholds for Enformer SAD track scores, and then adjust the p-value based upon the number of SNPs filtered from the given study. This is done with gwas_1_leading_snps_procedure.py for gwas 1, and gwas_3_scz_leading_snps_procedure.py for gwas 3, and gwas_4_alz_leading_snps_procedure.py for gwas_4.

Locations of Files


How to run specific files


Other files + Libraries



