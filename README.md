# bvr_zoops_code

Analysis and figure generation code for assessing spatial and temporal variability in zooplankton community structure

# Repo structure

This repository contains five main folders: 1) `conversions` contains files used to calculate zooplankton size from stage micrometer measurements, biomass using length-weight regressions (Downing and Rigler, 1984), and a taxonomic identification look-up table, 2) `figures` contains manuscript figures produced from various R scripts, 3) `output` contains files and tables summarizing statistical analyses, 4) `scripts` contains all code necessary to reproduce figures and analyses to assess spatial and temporal variability in zooplankton communities in Beaverdam Reservoir, and 5) `zoop_data` contains the raw zooplankton count and size data that are used to generate summary stats.

# Instructions to reproduce figures and analyses:

1.  Run `01_install.R` to install and load packages and functions necessary for R scripts.
2.  Run `02_VA_reservoirs_SummaryStats.R` in the `scripts` folder to generate summary statistics (i.e., density and biomass) for each zooplankton sample.
3.  Run `03_DVM_calcs.R` in the `scripts` folder to calculate hypolimnetic density and biomass from full and epilimnetic zooplankton tows.
4.  Run `04_All_MSN_zoop_figs.R` in the `scripts` folder to create the zooplankton normalized density figure in the manuscript (Fig. 3) and some additional SI density figures (Figs. S6, S11).
5.  Run `05_Migration_metrics.R` in the `scripts` folder to calculate metrics of DVM and DHM and generate Fig. 8, Fig. S9.
6.  Run `06_zoops_multivariate_stats.R` in the `scripts` folder to perform NMDS on zooplankton data, calculate metrics of variability (dispersion and pairwise), and create NMDS plots (Figs. 2, 4-7, Figs. S1-5, S7-8, S10, S12-14).

Note that `zoops_EDI.R` script was used to get zooplankton data into the correct format for publication on EDI (<https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=197&revision=3>).
