# bvr_zoops_code

Analysis and figure generation code for assessing spatial and temporal variability in zooplankton community structure

# Instructions to reproduce figures and analyses:

1.  Run `install.R` to install and load packages and functions necessary for R scripts.
2.  Run `VA_reservoirs2019_SummaryStats.R`, `VA_reservoirs2020_SummaryStats.R`, and `VA_reservoirs2021_SummaryStats.R` in the `scripts/annual_calcs` folder to generate summary statistics (i.e., density and biomass) for each zooplankton sample.
3.  Run `DVM_calcs-Summer19.R`, `DVM_calcs-Summer20.R`, and `DVM_calcs-Summer21.R` in the `scripts/annual_calcs` folder to calculate hypolimnetic density and biomass from full and epilimnetic zooplankton tows.
4.  Run `All_MSN_zoop_figs.R` in the `scripts` folder to create the zooplankton normalized density figure in the manuscript (Fig. 4) and some additional SI density figs (Figs. S3-4).
5.  Run `zoops_multivariate_stats.R` in the `scripts` folder to perform NMDS on zooplankton data, calculate metrics of variability (dispersion and pairwise), and create NMDS plots (Figs. 5-8, Fig. S1, Figs. S5-7).
6.  Run `Migration_metrics.R` in the `scripts` folder to calculate metrics of DVM and DHM and generate Fig. 9.
7.  Run `2020-2021_NMDS.R` to create the hourly driver multipanel figure using 2020-2021 zooplankton data (Fig. S8).

Note that `zoops_EDI.R` script was used to get zooplankton data into the correct format for publication on EDI (LINK HERE).
