# Aube_et_al_2022
Scripts for analyses, simulations and figure generation from Aubé et al. (2022).

Data preparation and figure generation steps are provided as Google Colaboratory notebooks, while Python scripts are included for the simulations and the early processing of the resulting data. 

## Data processing
The dataset from Hausser et al. (2019) is processed in Data_preparation.ipynb, in Data_preparation, to add P1-P2 associations as well as WGD/SSD annotations. The resulting files rates_Hausser.csv and couples_divergence.csv, in Data_sim_ready, were used to initialize all simulations. The Data_init subdirectory contains all the data files from external sources which are used in Data_preparation.ipynb.

## Simulations
All runs of *in silico* evolution were performed using the Genome_script.py script. This simulation program implements a few features which were not fully (or not all) explored in the paper, including loss-of-function mutations and gene-specific distributions of mutational effects (according to current expression levels).
The command lines corresponding to each of the calls of this simulation script are provided in commands_evol.txt. Comments have been added to separate the different groups of simulations.

The summary statistics for all simulations of the same group were compiled using the compile_data.py script.
Each of the calls of this script are provided in compiling_commands.txt.

The file evol_funct.py contains custom Python function used in Genome_script.py.

For each series of simulations, all data used in the construction of figures have been deposited:

* Data_full_dists: Full distributions of transcription and translation rates (for a small subset of simulations) shown in Fig S14.
* Full_data_WGD: Full evolutionary trajectories (transcription and translation rates through all mutation-selection rounds) for the small-scale 'mock' simulation shown in Fig 3.
* Sims_Bp_higher: Summary statistics and final log2-fold changes shown in Fig S8 for simulations assuming a larger mutational target size for translation rate.
* Sims_dupli_effect: Divergence correlations shown in Fig S9, from simulations investigating the effect of assumptions of the optimality of the immediate post-duplication expression levels.
* Sims_SSD: Summary statistics for simulations where the end condition was assessed according to yeast SSD-derived paralogs, shown in Fig S13.
* Summ_stats_constraints_simulations: Summary statistics and divergence correlations for simulations performed using asymmetrical or bivariate (with correlated transcription-translation effects) distributions of mutational effects. Results shown in Figs 5, S11 and S12.
* Summ_stats_dists: Summary statistics and divergence correlations resulting from simulations using normally distributed mutational effects (shown in Figs 4, S6, S7 and S15), as well as summary statistics for the screening of various standard deviations of mutational effects when using a bivariate distribution (Fig S10).  

## Statistical analyses and figure generation
Each figure shown in the article has been generated in a separate notebook, except for Figs 1 and S1 to S4 (combined in Fig1_Supps.ipynb). In each case, the notebook also includes all necessary statistical analyses and data processing (except for the calculation of summary statistics performed by the simulation script and their compilation using compile_data.py). All notebooks have been added in the Figures directory. Fig 2 is not included, as it was made in Inkscape. Similarly, Figs 4A and 5A were made using BioRender.com and are not included.

Any external data used to generate figures is included in the Analysis_v2 subdirectory.

## Versions
Simulations were performed using Python (v 3.10). All included stastistical tests and optimization routines were done using their SciPy (v 1.8.0) implementation.

The notebooks were executed in Google Colaboratory environments, with Python (v. 3.8.16) and SciPy (1.7.3).
