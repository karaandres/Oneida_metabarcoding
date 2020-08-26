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

**Analysis:**

    1. Quality control (MultiQC)
        - Script: multiqc.sh
        - Output files: multiqc_report.html

    2. Trim heterogeneity spacers: Trimming_HS_Metagenomics
        - Scripts: Trimming_HS_Primers.pl, forward reads: Running_Script_R1.sh, reverse reads: Running_Script_R2.sh

    3. Trim adaptors: Trimmomatic
       - Script: 
    
    4. Denoising and error removal (DADA2)
        - Script: dada2_peru_metabarcoding.sh
        - Output files:

    5. Assign taxonomic information to ASVs (BLASTn)
        - Used tax_trace.pl to obtain higher taxonomic information from staxids
        - Scripts:
        - Output files:

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
