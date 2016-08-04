R scripts for analysis of Musca immunity
==

These are the R scripts and raw data associated with the manuscript Sackton et al ()....

To replicate our analysis:

1. Load data and run the differential expression analysis with DESeq2 by running the load_data.R script (which calls difexp.R). Note that this script has a number of paths that will need to be adjusted to run properly.

2. Run analysis scripts:

Results section 1, "Identifying genes regulated by infection in Musca domestica":
	difexp_analysis.R

Results section 2, "Gene ontology analysis suggests a coordinate shift...":
	GOanalysis.R

Results section 3, "Comparison to D. melanogaster RNA-seq data...":
	difexp_dmel_vs_musca.R

Results, section 4, "The infection-induced transcriptome of Musca domestica is enriched...":
	strata_analysis.R

3. Prepare duplication and loss data for analysis by running the clean_duploss_data.R script.

Results, section 6, "Poisson regression is an accurate method for estimating rates of gene gain and loss":

a. Estimating rates from simulated data: duploss_simulation_analysis.R
b. Testing Poisson mixed model for well-calibrated false positive rates: 
	In the R/modeltest directory, first generate simulation results with:
	--compute_falsepos_bysim.R
	--compute_falseneg_bysim.R
	
	Then analyze results with:
	--analyze_sim_results.R
	
	Note that these scripts run some other models that are not discussed in the manuscript, as the best approach is the mixed model approach.
	
4. Analysis of duplication and loss results.

Results, section 7, "Genes encoding effector and recognition proteins duplicate rapidly...": 
	duploss_analysis.R
	
5. Generate figures and tables for paper:
	- figures.R generates all the figures for our manuscript
	- tables.R generates data for tables that are not already generated in other scripts