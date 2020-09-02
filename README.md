# Oneida Lake: sampling fishes using eDNA metabarcoding vs. traditional gears

## This repository contains the files related to the data processing and data analysis for eDNA metabarcoding of fishes in Oneida Lake, NY

### Overview

In this project, we test the effectiveness of using environmental DNA (eDNA) approaches to characterize a fish community compared to historical fish records and traditional sampling gears including electrofishing, seining, and  We target a hypervariable region of the 12S rRNA gene to characterize fish communities across 25 sites in Oneida Lake.

### Methodological approach

**eDNA sampling:**
  - At each sampling site, we collected 2 L surface water samples
  - Environmental measurements characterizing the hydrology and water chemistry of the rivers were taken at each sample point
  - DNA was extracted and PCR-amplified using a fish-specific mitochondrial 12S rDNA primer pair (MiFish, Miya et al.).

**Sequencing:**
  - Illumina NextSeq 300 bp Mid output kit
  - Raw read output: 124763444

**Analysis:**
  1. Trim heterogeneity spacers: [Trimming_HS_Metagenomics](https://github.com/noushing/Trimming_HS_Metagenomics)
     - Scripts: [Trimming_HS_Primers.pl](Trimming_HS_Primers.pl), forward reads: [Running_Script_R1.sh](Running_Script_R1), reverse reads: [Running_Script_R2.sh](Running_Script_R2.sh)

  2. Trim adaptors: Trimmomatic
     - Script: [trimmomatic_oneida.sh](trimmomatic_oneida.sh)
     - Reads remaining after trimming spacers and adaptors: 94717624
     
  3. Quality control (MultiQC)
     - Script: [multiqc.sh](multiqc.sh)
    
  4. Denoising and error removal (DADA2)
     - Script: [dada2_oneida_metabarcoding.txt](dada2_oneida_metabarcoding.txt)
     - Output files: [seqtab.nochim.csv](seqtab.nochim.csv)

  5. Assign taxonomic information to ASVs (BLASTn)
     - Used BLAST and [tax_trace.pl](https://github.com/theo-allnutt-bioinformatics/scripts/blob/master/tax_trace.pl) to obtain higher taxonomic information for top 5 BLAST hits
     - Assigned species to ASV when top species hit was > 98% identity. Ambiguous hits or hits not exceedning 98% identity were assigned higher taxonomic information. 
     - Scripts: 
       - [blastn_oneida.txt](blastn_oneida.txt)
       - [assigning_higher_taxonomic_info.R](assigning_higher_taxonomic_info.R)
       - [ASV_to_spp.R](ASV_to_spp.R)
     - Output files: 
       - [blstn_nonchim_fmt10_scinames_Oneida_withseq_taxonimic_matched.csv](blstn_nonchim_fmt10_scinames_Oneida_withseq_taxonimic_matched.csv)
       - [ASV_count_by_sites.csv](ASV_count_by_sites.csv)
       - [sp_count_by_site.csv](sp_count_by_site.csv)

  6. Data filtering
     - Removed ASVs with low read counts and non-target taxa
     - Scripts:
     - Output files:

  7. Diversity metrics and site comparisons
     - Scripts:
     - Output files:

  8. Comparison to species list
     - Compare species matches to traditional sampling techniques
     - Scripts:
     - Output files:
