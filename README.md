# Aube_et_al_2022
Scripts for analyses, simulations and figure generation from Aube et al. (2022).

Data preparation and figure generation steps are provided as Google Colaboratory notebooks, while Python scripts are included for the simulations and the early processing of the resulting data. 

## Data processing
The dataset from Hausser et al. (2019) is processed in Data_preparation.ipynb, in Data_preparation, to add P1-P2 associations as well as WGD/SSD annotations. The resulting files rates_Hausser.csv and couples_divergence.csv, in Data_sim_ready, were used to initialize all simulations. 

## Simulations
All runs of *in silico* evolution were performed using the Genome_script.py script. This simulation program implements a few features which were not fully (or not all) explored in the paper, including loss-of-function mutations and gene-specific distributions of mutational effects (according to current expression levels). 

The file evol_funct.py contains custom Python function used in Genome_script.py.

For each series of simulations, the corresponding batch files as well as any Python script used to pre-process the simulation results were combined into one directory. In addition, all files used for the generation of the figures, containing raw or processed data, have been included.

* Mock_WGD = 'Mock' simulation of 50 paralog pairs shown in Fig 3. **Data not included (files too big)**
* Mut_SD_WGD = Series of simulations under the standard framework to identify the best-fitting standard deviation of mutational effects under a high efficacy of selection (N = 1e6). Data used for Figs 4, S3 and S4. **Data not included (files too big)**
* Mut_SD_WGD_1e5 = Same as above, but for N = 1e5. Data used for Figs 4 and S7. **Data not included (files too big)**
* Mut_alpha_WGD = Simulations using a skew normal distribution of mutational effects, with various values of the skewness parameter alpha. Data used in Figs 5 and S8.
* Mut_alpha_WGD_1e5 = Same as above, but for N = 1e5. Data used in Fig S9.
* Mut_SD_bivariate_WGD = Series of simulations performed using a bivariate distribution of mutational effects to identify the best-fitting "reference" standard deviation of effects, under a high efficacy of selection (N = 1e6). Data used for Fig S7.
* Mut_SD_bivariate_WGD_1e5 = Same as above, but for N = 1e5. Data used in Fig S7.
* Mut_bivariate_WGD = Series of simulations to test various levels of correlation between transcriptional and translational mutations when a bivariate distribution of mutational effects is used. Performed under N = 1e6 and using the corresponding best-fitting "reference" standard deviation of mutational effects. Data used in Figs 5 and S8.
* Mut_bivariate_WGD_1e5 = Same as above, but for N = 1e5 (and the corresponding best-fitting "reference" standard deviation of mutational effects. Data used in Fig S9.
* Bp_higher_WGD = Series of simulations allowing a larger mutational target size for translation rate. Performed using the standard framework and under an assumption of high selection efficacy (N = 1e6). Data used for Fig S5.
* Dupli_effect_WGD = Series of simulations performed while varying the post-duplication change of the optimum of cumulative protein abundance, using the standard framework under N = 1e6. Data used for Fig S6.
* Mut_SD_SSD = Series of simulations under the standard framework to identify the best-fitting standard deviation of mutational effects under a high efficacy of selection (N = 1e6), but this time when trying to replicate the divergence patterns of SSD-derived paralogs. Data used in Fig S10.
* Mut_SD_SSD_1e5 = Same as above, but for N = 1e5. Also used for Fig S10.    

## Statistical analyses and figure generation
Each figure shown in the article has been generated in a separate notebook, except for Figs 1 and S1 (combined in Fig1_Supp1.ipynb). In each case, the notebook also includes all necessary statistical analyses and data processing (except the pre-processing done using Python scripts included in the simulation directories). All notebooks have been added in the Figures directory. Fig 2 is not included, as it was made in Inkscape.
