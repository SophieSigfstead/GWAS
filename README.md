# GWAS Repository Information
## Overall Protocol

This document provides some details on the repository. This repository contains several manual trials of our GWAS SNP filtering method. 

The main protocol for each GWAS analysis is to:
1. Identify an Enformer track list relevant to the trait being studied. Tracks are labelled according to tissue type in target_dnase_ataq_tracks_labelled.csv.

2. Merge 1000 genomes SAD scores for each track in the track list with summary statistics for the given GWAS using a gwas_#_matching.py. These files match rows of summary statistics by either snp id (rs id) or chromosome and position pairs, depending on the information provided in the summary statistics. For gwas 1, this is called "gwas_1_matching.py". For gwas 3, this is called "gwas_3_scz_matching.py". For gwas 4, this is called "gwas_4_matching.py". 

To replicate this in a new study: 
The 1000 genomes SAD score data are stored in 1000genomes_as_csv. The summary statistics for your given study must be downloaded to the repository. Then, the SAD score data must be merged with the summary statistics on a common column depending on the variables in the summary stats (e.g. rs-id, [chr,pos]). 

3. Replicate the given GWAS using a pre-filtered list of SNPs. These files are labelled gwas_#_leading_snps_procedure.py. These files first filter the SNPs based upon global thresholds for Enformer SAD track scores, and then adjust the p-value based upon the number of SNPs filtered from the given study. This is done with gwas_1_leading_snps_procedure.py for gwas 1, and gwas_3_scz_leading_snps_procedure.py for gwas 3, and gwas_4_alz_leading_snps_procedure.py for gwas 4.

To replicate this in a new study, 1st apply GWAS filtering. This is done by computing global threshhold for Enformer SAD track scores across n_0 initial SNPs, then filtering the SNP list according to these thresholds. Compute the new number of SNPs (n_1) and multiply 5e-8 x (n_0/n_1). This will create a higher p-value threshold. Then recreate the GWAS selection process using methods outlined in the study using this new threshold for significance and filtered snp list. 

4. Following the above steps, perform analyses of results compared to original study. 

Note that all files associated with a specific gwas are labelled beginning with “gwas_#”. To reference which gwas study this corresponds to please look at references.txt. Note that the analysis for gwas 2 was not completed due to limitations in study replicability. The following analyses have been conducted, as of Oct 30, 2024:

### GWAS 1:

- Original summary statistics are available in GWAS/gwas_1_and_2_summary_statistics_data/GCST90277450.tsv
- The matched files (result of step 2 in Overall Protocol) are stored in gwas_1_matching. The file for a track t is labelled as result_SAD<t>.csv, and houses the SAD scores for each SNP in 1KG. This can be reproduced by running: python gwas_1_matching.py ./gwas_1_and_2_summary_statistics_data/GCST90277450.tsv ./1000genomes_as_csv [0,1,9,76,78,80,81,172,179,216,240,261,278,319,326,338,355,370,403,411,421,458,462,469,499,524,552,580,582,602,644,669]
- Files produced as gwas_1_leading_SNPs_procedure.py runs are stored in GWAS/gwas_1_intermediate_files
- Results of running gwas_1_leading_snps_procedure.py of each threshold t are stored in GWAS/filtered_snps_gwas_1/filtered_snps_gwas_1_threshold=<t>.csv.

Analysis of results of GWAS 1 analysis was completed using the following data:
- eQTL (Data available in GWAS/EQTLs)
- FANTOM5_annotations (Data available in GWAS/FANTOM5_annotations)
- PAINTOR

- A seperate analysis was ran to understand the effects of a single SAD track being used as a threshold. The results of this are stored in GWAS/gwas_1_single_track_analysis/GWAS_1_leading_SNPs_by_track amd GWAS/gwas_1_single_track_analysis/GWAS_1_leading_snps_by_track_random. These results were produced using the GWAS/gwas_1_single_track_analysis/gwas_1_leading_SNPs_by_track.py file. The comparison notebook to all results is labelled GWAS/gwas_1_single_track_analysis/gwas_1_all_vs_single_track.ipynb and for  a random set of tracks in GWAS/gwas_1_single_track_analysis/gwas_1_all_vs_single_track_random.ipynb. A large-scale version of this analysis is currently underway. 

### GWAS 2:
- Original summary statistics are available in GWAS/gwas_1_and_2_summary_statistics_data/PGC_UKB_depression_genome-wide.txt
- - The matched files (result of step 2 in Overall Protocol) are stored in gwas_1_matching. The file for a track t is labelled as result_SAD<t>.csv, and houses the SAD scores for each SNP in 1KG. This can be reproduced by running: python gwas_2_matching.py ./gwas_1_and_2_summary_statistics_data/PGC_UKB_depression_genome-wide.txt ./1000genomes_as_csv [0,1,9,76,78,80,81,172,179,216,240,261,278,319,326,338,355,370,403,411,421,458,462,469,499,524,552,580,582,602,644,669]
- Results of gwas_2_leading_SNPs_procedure_v2.py for each threshold t are stored in GWAS/filtered_snps_gwas_2/filtered_snps_gwas_2_threshold=<t>.csv.
- Files produced as gwas_2_leading_SNPs_procedure_v2.py runs are stored in GWAS/gwas_2_intermediate_files.

### GWAS 3:
- Summary statistics are available in GWAS/gwas_3_scz_original_files
- Files produced as gwas_3_scz_leading_snps_procedure.py runs are stored in GWAS/gwas_3_scz_intermediate_files.
- Results of running gwas_3_scz_leading_snps_procedure.py are available in GWAS/gwas_3_scz_result_files. 

### GWAS 4:
- Summary statistics are available in GWAS/gwas_4_alz_summary_statistics. The study we are looking at is labelled PGCALZ2sumstatsExcluding23andMe.txt. The stats from a 2019 ALZ GWAS study from the same group (precursor to this one) are also there as they are used to label snps in this study with the correct snp ids. 
- The original results of the study, which included 23andMe data (we don't have access to) are available in gwas_4_alz_original_results_final.csv
- The matched files (result of step 2 in Overall Protocol) are stored in gwas_1_matching. The file for a track t is labelled as result_SAD<t>.csv, and houses the SAD scores for each SNP in 1KG. This can be reproduced by running: python gwas_4_alz_matching.py ./GWAS/gwas_4_alz_summary_statistics/PGCALZ2sumstatsExcluding23andMe.txt ./1000genomes_as_csv [0,1,9,76,78,80,81,172,179,216,240,261,278,319,326,338,355,370,403,411,421,458,462,469,499,524,552,580,582,602,644,669].
- Files produced as gwas_4_alz_leading_snps_procedure.py runs are stored in GWAS/gwas_4_alz_intermediate_files.
- Results of running gwas_4_alz_leading_snps_procedure.py are available in GWAS/gwas_4_alz_result_files

## Other Directories:
### GWAS/enformer_and_tangermeme 
- At beginning of project, wanted to ensure that we could run Enformer + Tangermeme in pytorch. Files to do so are contained here. These were just tests to ensure functionality, nothing is built yet. 

### Finemapping: GWAS/finemapping and GWAS/fine_mapping_results
- Used to do finemapping analysis on GWAS 1

### Genome Wide Complex Trait Analysis (GCTA): GWAS/gcta_output_directory and GWAS/gcta
- GCTA software website: https://yanglab.westlake.edu.cn/software/gcta/#Overview
- GCTA install is in GWAS/gcta
- GWAS/gcta_output_directory housed intermediate gcta files while gwas_2 was running, but is currently empty. GWAS 2 was not able to be replicated well, so the analysis was not completed. 

### gencode_annotations
- Contains annotation from: https://www.gencodegenes.org/human/release_19.html
- Used in some analysis in the project (?), however the annotation in genome_assembly was more commonly used and up to date

### genome_assembly
- Contains genome annotations. Version v46lift37 (GWAS/genome_assembly/gencode.v46lift37.annotation.gff3) is most applicable to hg37, which is applicable to gwas 1-4.
- Contains make_coding_regions_list.py. This was used to create the exon_regions that are used to identify coding snps in the GWAS filtering procedure. This routine creates the regions for coding snps, which are contained in GWAS/genome_assembly/exon_regions_v2.csv

### h5_to_csv_scripts
- Used to convert h5 files to csv. This was used to produce the 1000genomes_as_csv file from H5 files provided by Enformer paper: 
https://console.cloud.google.com/storage/browser/dm-enformer/variant-scores/1000-genomes/enformer;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

### leading SNPs_SAD_track_visualization
- This was an analysis to compare the SAD enformer scores of the track list we identified for psychiatric conditions ([0,1,9,76,78,80,81,172,179,216,240,261,278,319,326,338,355,370,403,411,421,458,462,469,499,524,552,580,582,602,644,669]) to a random track list (see gwas_random_tracks_matching for the track list and files analysed).

### random_intermediate_files
- Used for running analyses for Gwas 1 on random track sets, any intermediate files were placed here instead of gwas_1_intermediate_files



