# bvr_zoops_code

Analysis and figure generation code for assessing spatial and temporal variability in zooplankton community structure

# Instructions to reproduce figures and analyses:

1.  Download data from EDI (*link here once published!*)
2.  Run `install.R` to install and load packages and functions necessary for R scripts.
3.  Run `VA_reservoirs_SummaryStats.R` in the `scripts/annual_calcs` folder to generate summary statistics (i.e., density and biomass) for each zooplankton sample.
4.  Run `DVM_calcs.R` in the `scripts/annual_calcs` folder to calculate hypolimnetic density and biomass from full and epilimnetic zooplankton tows.
5.  Run `All_MSN_zoop_figs.R` in the `scripts` folder to create the zooplankton normalized density figure in the manuscript (Fig. 4) and some additional SI density figs (Figs. S3-4).
6.  Run `zoops_multivariate_stats.R` in the `scripts` folder to perform NMDS on zooplankton data, calculate metrics of variability (dispersion and pairwise), and create NMDS plots (Figs. 5-8, Fig. S1, Figs. S5-7).
7.  Run `Migration_metrics.R` in the `scripts` folder to calculate metrics of DVM and DHM and generate Fig. 9.
8.  Run `2020-2021_NMDS.R` to create the hourly driver multipanel figure using 2020-2021 zooplankton data (Fig. S8).

Note that `zoops_EDI.R` script was used to get zooplankton data into the correct format for publication on EDI (LINK HERE).

## Data availability

Zooplankton data used in this analysis are published to the EDI data repository (link here). Coversions are in the conversions folder, etc.