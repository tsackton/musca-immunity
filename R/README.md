R scripts for analysis of Musca immunity
==

These are the R scripts and raw data associated with the manuscript Sackton et al ()....

To replicate our analysis:

1. Load data and run the differential expression analysis with DESeq2 by running the load_data.R script (which calls difexp.R)

2. Run analysis scripts:
  - GOanalysis.R runs the analysis of GO terms, using blast2go data for M. domestica from the genome paper.
  - strata_analysis.R runs the phylostratigraphic analysis.
  - difexp_analysis.R to run the basic analysis of differential expression classes
  - analyze_new_models.R to calculate some basic stats on the new PASA gene models

3. Run gene family evolution scripts:
	- clean_duploss_data.R creates the ultrametric tree and prepares data for both CAFE and the Poisson modeling described in the manuscript.
	- per_ogs_rates.R estimates rates for each OGS using the Poisson model and then creates rates for simulation
	- duploss_simulation_analysis.R runs the analysis of the simulated gene families
	- duploss_analysis.R runs the main code to analyze duplication and loss rates, as well as export additional data for CAFE

4. Generate figures for paper:
	- figures.R generates all the figures for our manuscript