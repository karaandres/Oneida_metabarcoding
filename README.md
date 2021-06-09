# Oneida Lake: sampling fishes using eDNA metabarcoding vs. traditional gears

## This repository contains the files related to the data processing and data analysis for eDNA metabarcoding of fishes in Oneida Lake, NY

### Overview

In this project, we test the effectiveness of using environmental DNA (eDNA) approaches to characterize a fish community compared to historical fish records and traditional sampling gears including electrofishing, seine nets, fyke nets, and gill nets. We target a hypervariable region of the 12S rRNA gene to characterize fish communities across 25 sites in Oneida Lake.

### Methodological approach

**eDNA sampling:**
  - At each sampling site, we collected 2 L surface water samples
  - Environmental measurements characterizing the hydrology and water chemistry of the rivers were taken at each sample point
  - DNA was extracted and PCR-amplified using a fish-specific mitochondrial 12S rDNA primer pair (MiFish, Miya et al. 2015).

**Sequencing:**
  - Illumina NextSeq 300 bp Mid output kit
  - Raw read output: 124763444

**Bioinformatic processing:**
  1. Trim heterogeneity spacers: [Trimming_HS_Metagenomics](https://github.com/noushing/Trimming_HS_Metagenomics)
     - Script: [Trimming_HS.txt](scripts/Trimming_HS.txt)
     - Reads remaining: 103024059

  2. Trim adaptors and reads < 35 bp: Trimmomatic
     - Script: [trimmomatic_oneida.sh](scripts/trimmomatic_oneida.sh)
     - Reads remaining: 94717624
     
  3. Quality control (MultiQC)
     - Script: [MultiQC.txt](scripts/MultiQC.txt)
    
  4. Denoising and error removal (DADA2)
     - Script: [dada2_oneida_metabarcoding.txt](scripts/dada2_oneida_metabarcoding.txt)
     - Output file: seqtab.nochim.csv (total 2996 ASVs)

  5. Assign taxonomic information to ASVs (BLASTn)
     - Used BLAST and [tax_trace.pl](https://github.com/theo-allnutt-bioinformatics/scripts/blob/master/tax_trace.pl) to obtain higher taxonomic information for top 5 BLAST hits
       - [blastn_oneida.txt](scripts/blastn_oneida.txt)
     - Selected taxonomic information of interest and rearranged in workable format 
       - Script: [assigning_higher_taxonomic_info.R](scripts/assigning_higher_taxonomic_info.R)
       - Output file: blstn_nonchim_fmt10_scinames_Oneida_withseq_taxonimic_matched.csv
     - Assigned species to ASV when top species hit is >= 98% identity. Ambiguous hits or hits not exceeding 98% identity were assigned higher taxonomic information. Taxonomic assignments were added to the matrix of read counts of ASVs per site.
       - Script: [ASV_to_spp_oneida.R](scripts/ASV_to_spp_oneida.R)
       - Output file: [ASV_count_by_sites.csv](datasets/ASV_count_by_sites.csv)

  6. Data filtering
     - Removed ASVs with low read counts (<0.1% of all reads in a sample), removed non-target taxa, and adjusted reads to account for false positives in control samples (removed ASVs with read counts < the mean number of non-zero reads across all negative controls)
     - Removed negative control samples and collapsed ASV read counts per site by species
       - Script: [data_filtering_oneida.R](scripts/data_filtering_oneida.R)
       - Output file: [sp_read_count_by_site.csv](datasets/sp_read_count_by_site.csv)

  7. Combine eDNA and traditional gear datasets
     - Standardized all fish survey datasets (eDNA, EF, seine, fyke, gillnet, and historical)
     - Outputs: 1) all gears separate files: [eDNA_dat.csv](datasets/eDNA_dat.csv), [ef_dat.csv](datasets/ef_dat.csv), [fyke_dat.csv](datasets/fyke_dat.csv), [gillnet_dat.csv](datasets/gillnet_dat.csv), [seine_dat.csv](datasets/seine_dat.csv); 2) all gears combined files: presence-absenece ([all_methods_presence.csv](datasets/all_methods_presence.csv)), abundance ([all_methods_abundance.csv](datasets/all_methods_abundance.csv))

**Data analysis:**
  1. Presence-absence comparisons at 3 spatial scales: lake-wide, habitat, and site 
     - Contains: lake-wide species overlap among all datasets, lake-wide species accumulation curves and richness estimates, comparisons of species occurrence for eDNA vs. traditional gears across pelagic and nearshore habitats, NMDS analyses for eDNA and traditional gears across habitat types, site-level comparisons of species richness for eDNA vs. electrofishing
     - Script: [Oneida_presence_absence_comparisons.R](scripts/Oneida_presence_absence_comparisons.R)
     
  3. Species abundance comparisons at 3 spatial scales: lake-wide, habitat, and site 
     - Script: 
