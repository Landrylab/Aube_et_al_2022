#!/usr/bin/env python
import argparse
import datetime
import time
import os
import evol_funct  # Custom functions
import scipy.stats as stats
import pandas as pd
import numpy as np
import scipy.optimize as optimize
import math
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
import statannot

PARSER = argparse.ArgumentParser(description="Simulating the divergence of duplicate genes following a duplication "
                                             "event under quantitative subfunctionalization while optimizing the "
                                             "tradeoff between the economy and the precision of gene expression.")

PARSER.add_argument("-n", "--name",
                    metavar="STRING", dest="run_name", default=f'{str(datetime.date.today())}', type=str,
                    help="Prefix under which the current run will be saved (defaults to current date, in "
                         "YYYY-MM-DD format).")

PARSER.add_argument("--path_out",
                    metavar="STRING", dest="path_out", default='', type=str,
                    help="Path to the directory where the 'set' and 'folder' folders should be created (or "
                         "are saved, if they already exist).")

PARSER.add_argument("--set",
                    metavar="STRING", dest="sim_set", default='test_set', type=str,
                    help="Identifier for the larger set of simulations of which the current run is part. "
                         "Will also be the name of the main folder in which the subfolder is created.")

PARSER.add_argument("--folder",
                    metavar="STRING", dest="sim_folder", default='test_folder', type=str,
                    help="Identifier of the subfolder (within the main folder defined above) in which all output "
                         "files will be saved.")

PARSER.add_argument("--path_data",
                    metavar="STRING", dest="path_data", default='', type=str,
                    help="Path to the files rates_Hausser.csv and couples_divergence.csv containing input data.")

PARSER.add_argument("--make_figures",
                    metavar="BOOL", dest='make_figs', default=False, type=bool,
                    help="Boolean specifying whether the figures should be created and saved at the end of the "
                         "simulation run. Defaults to False.")

PARSER.add_argument("-c", "--couples",
                    metavar="INT", dest="n_couples", default=2000, type=int,
                    help="Number of paralog couples simulated (defaults to 2000).")

PARSER.add_argument("--seed",
                    metavar="INT", dest="seed", default=None, type=int,
                    help="Seed used to initialize the random number generator.")

PARSER.add_argument("-o", "--loss_of_function",
                    metavar="BOOL", dest="lof", default=False, type=bool,
                    help="Specifies whether gene loss is allowed during the evolution of the duplicate couples. When"
                         " set to True, the loss of a copy happens as soon as it becomes effectively neutral. Defaults"
                         " to False.")

PARSER.add_argument("--pop",
                    metavar="FLOAT", dest="effective_pop", default=1000000, type=float,
                    help="Effective population size to be considered for the computation of mutation fixation "
                         "probabilities (defaults to 1000000). Scientific notation (ex: 1.0e6) can be used as well.")

PARSER.add_argument("--sampling",
                    metavar="{random, genes}", dest="samp_meth", default="random", type=str,
                    choices=["random", "genes"],
                    help="Method used to sample ancestral gene to initialize the simulation. If set to 'random', all "
                         "parameters are drawn for the corresponding distributions (or set to an arbitrary value), "
                         "while if set to 'genes', S. cerevisiae genes are sampled at random to preserve relationships "
                         "between expression rates (Bm and Bp) and Q as well as lm. Defaults to 'random'.")

PARSER.add_argument("--q_value",
                    metavar="FLOAT", dest="q_val", default=None, type=float,
                    help="Only relevant when ancestral genes sampling method is set to 'random'. Allows to set "
                         "an arbitrary value for the noise sensitivity (curvature of the fitness function) for "
                         "all replicate ancestral genes. Defaults to None, in which case the complete distribution of "
                         "Q values inferred from Hausser et al. (2019) will be used.")

PARSER.add_argument("--lm_value",
                    metavar="FLOAT", dest="lm_val", default=1350, type=float,
                    help="Only relevant when ancestral genes sampling method is set to 'random'. Allows to set an "
                         "arbitrary value for the pre-mRNA lengths of all replicate ancestral genes. Defaults to 1350, "
                         "the median value reported by Hausser et al. (2019). If set to None, the complete "
                         "distribution of lengths from Hausser et al. (2019) is used.")

PARSER.add_argument("--end",
                    metavar="{pEst, loss}", dest="end_cond", default="pEst", type=str,
                    choices=["pEst", "loss"],
                    help="End condition of the simulation. Specifies whether is it stopped when a realistic level of"
                         " protein abundance divergence is reached or when a sufficient (realistic) number of copies "
                         "have been lost. Default to 'pEst'.")

PARSER.add_argument("-v", "--cv0",
                    metavar="FLOAT", dest="cv_0", default=0.1, type=float,
                    help="Value of the noise floor (cv0), which represents the extrinsic component of protein "
                         "expression noise. Defaults to 0.1, the value reported by Hausser et al. (2019) for yeast.")

PARSER.add_argument("-u", "--dupli_effect",
                    metavar="FLOAT", dest="dupli_effect", default=1.87, type=float,
                    help="Variation of the optimal protein abundance after the duplication. Defaults to 1.87, meaning"
                         " that the immediate post-duplication expression level is not optimal. Should be set between "
                         "1.87 and 2.14 (approximate values) to ensure positive post-duplication fitness for all "
                         "simulated pairs.")

PARSER.add_argument("--dupli_m",
                    metavar="FLOAT", dest="total_bm", default=2.0, type=float,
                    help="Float representing the post-duplication variation of the total transcription rate. Defaults "
                         "to 2.0, meaning that both duplicates have the ancestral transcription rate and that the"
                         " total transcription is thus doubled.")

PARSER.add_argument("-m", "--mut_ratio",
                    metavar="FLOAT", dest="mut_ratio", default=1, type=float,
                    help="Ratio of mutational target sizes for transcription and translation rates "
                         "(transcription / translation). Defaults to 1, meaning equal mutation probability for "
                         "both traits")

PARSER.add_argument("--bm_dist",
                    metavar="(mean, std)", dest="bm_params", default='(0, 0.025)', type=str,
                    help="Parameters (mean and standard deviation) of the normal distribution of mutational effects for"
                         " transcription rate. Both parameters are expressed as a percentage of current expression."
                         " Defaults to '(0, 0.025)'.")

PARSER.add_argument("--bp_dist",
                    metavar="(mean, std)", dest="bp_params", default='(0, 0.025)', type=str,
                    help="Parameters (mean and standard deviation) of the normal distribution of mutational effects for"
                         " translation rate. Both parameters are expressed as a percentage of current expression."
                         " Defaults to (0, 0.025).")

PARSER.add_argument("--skew",
                    metavar="FLOAT", dest="skew", default=0, type=float,
                    help="Alpha skewness parameter for the skew normal distributions of mutational effects. Defaults to"
                         "0, meaning that both distributions are normal. This alpha parameter is also applied to both"
                         "variables in the bivariate implementation of the (skew normal) mutational effect "
                         "distributions.")

PARSER.add_argument("--bivariate",
                    metavar="BOOL", dest='bivariate', default=False, type=bool,
                    help="Boolean specifying whether transcriptional and translational mutational effects should be"
                         "drawn from a bivariate distribution. Defaults to 'False', meaning that a given mutation only"
                         "affects one of the two parameters and that mutational effects are drawn from two separate"
                         "distributions.")

PARSER.add_argument("--correlation",
                    metavar="FLOAT", dest='corr', default=0, type=float,
                    help="Correlation coefficient (Pearson) between the transcriptional and translational effects"
                         "of mutations. Defaults to 0. Only used if mutational effects are sampled from a bivariate"
                         "distribution.")

PARSER.add_argument("--bidirectional",
                    metavar="BOOL", dest="bidirectional", default=False, type=bool,
                    help="Boolean specifying whether two mutational effects distribution with opposite skewness values"
                         "should be used for high and low values of transcription and translation rate. Defaults to "
                         "false. meaning that the alpha value specified is used for the sampling of all mutations. If"
                         "set to True, a positive skewness is used for expression rates below a certain threshold,"
                         "while a negative skewness is used for values above a threshold. Currently, this is only"
                         "implemented for mutations drawn from two univariate distributions.")

PARSER.add_argument("--threshold",
                    metavar="FLOAT", dest="ext_quant", default=0.1, type=float,
                    help="Quantile of expression rates value past which an asymmetrical mutational effects distribution"
                         "is used. A positive skewness is applied below this threshold quantile, while a negative"
                         "skewness is applied above the (1-threshold) quantile.")

PARSER.add_argument("-d", "--decay",
                    metavar="(alpha_m, alpha_p)", dest="decay_rates", default=(5.1, 1.34), type=tuple,
                    help="mRNA (alpha_m) and protein (alpha_p) degradation rates. Defaults to the typical decay rates"
                         " reported by Hausser et al. (2019), in h^-1.")

PARSER.add_argument("--bm_bounds",
                    metavar="(min, max)", dest="bm_bounds", default=None, type=str,
                    help="Tuple of minimal and maximal allowed transcription rates (in mRNAs per hour). If no bounds "
                         "are specified (default=None), they are respectively set to 0 and to the maximal value "
                         "reported by Hausser et al. (2019).")

PARSER.add_argument("--bp_bounds",
                    metavar="(min, max)", dest="bp_bounds", default=None, type=str,
                    help="Tuple of minimal and maximal allowed translation rates (in proteins per mRNA per hour). If no"
                         " bounds are specified (default=None), they are respectively set to 0 and to the maximal value"
                         " reported by Hausser et al. (2019).")

PARSER.add_argument("--dupli_type",
                    metavar="{WGD, SSD, all}", dest="dupli_type", default='all', type=str,
                    choices=["WGD", "SSD", "all"],
                    help="Argument to specify which type of duplicates should be used for the Mood's test of protein"
                         "abundance divergence, which is used for the 'pEst' end condition, Defaults to all, meaning"
                         "that all yeast duplicate couples will be used.")

PARSER.add_argument("--dataset",
                    metavar="{Original, Corrected}", dest="dataset", default="Original", type=str,
                    choices=["Original", "Corrected"],
                    help="Argument specifying whether the original or the 'corrected' translation rates from Hausser et"
                         "al. (2019) should be used to generate ancestral genes as well as for the comparisons.")

PARSER.add_argument("--epistasy",
                    metavar="{add, mult}", dest="epistasy", default="mult", type=str,
                    choices=["add", "mult"],
                    help="Argument specifying what kind of epistasy is assumed between mutations. Defaults to"
                         "multiplicative. If additive, the absolute magnitude of the mutational effects is set "
                         "according to the ancestral values of transcription and translation rates. If multiplicative,"
                         "the absolute value of mutational effects is set according to the current transcription or"
                         "translation rate at each mutation round.")

PARSER.add_argument("--ancestral",
                    metavar="{all, dupli, single}", dest="anc_dist", default="all", type=str,
                    choices=["all", "dupli", "single"],
                    help="Argument specifying whether the properties of the simulated ancestral singletons will be "
                         "drawn from the whole yeast genome, from the genes which are part of WGD- or"
                         "SSD-derived paralog couples or from the singletons. Defaults to 'all', meaning that the "
                         "complete genome is used.")

PARSER.add_argument("--full_data",
                    metavar="BOOL", dest="full_data", default=False, type=bool,
                    help="Boolean specifying whether the transcription and translation rates should be saved after "
                         "every mutation round. Defaults to False, meaning that only the initial and final rates"
                         "will be saved. Only set to True for runs made with a <500 duplicate couples.")

ARGS = PARSER.parse_args()
# Reference for total process time
t_start = time.process_time()

# 1- Arguments parsing and definition of cellular constants and gene properties distributions.
# A) Script arguments
run_name = ARGS.run_name
path_out = ARGS.path_out
sim_set = ARGS.sim_set
sim_folder = ARGS.sim_folder
path_data = ARGS.path_data
make_figs = ARGS.make_figs
full_data = ARGS.full_data

# B) The specified directories are created if they do not already exist
path_set = os.path.join(path_out, sim_set)
path_folder = os.path.join(path_set, sim_folder)

if not os.path.exists(path_set):
    try:
        os.makedirs(path_set)
    except FileExistsError:
        pass

if not os.path.exists(path_folder):
    try:
        os.makedirs(path_folder)
    except FileExistsError:
        pass

# C) Simulation arguments
n_couples = ARGS.n_couples
rng_seed = ARGS.seed
pop_size = ARGS.effective_pop
samp_meth = ARGS.samp_meth
q_val = ARGS.q_val
lm_val = ARGS.lm_val
dupli_loss = ARGS.lof
end_cond = ARGS.end_cond
dupli_effect = ARGS.dupli_effect
bm_bounds = ARGS.bm_bounds
bp_bounds = ARGS.bp_bounds
dupli_type = ARGS.dupli_type
dataset = ARGS.dataset
epistasy = ARGS.epistasy
anc_dist = ARGS.anc_dist
correlation = ARGS.corr
bidirectional = ARGS.bidirectional
ext_quant = ARGS.ext_quant

# D) Cellular constants

# Degradation rates
alpha_m = ARGS.decay_rates[0]
alpha_p = ARGS.decay_rates[1]

# Noise floor
cv_0 = ARGS.cv_0

# Transcription cost per nucleotide
c_m = 1.2e-09

# E) Gene properties distributions
# Initialization of the random number generator that will be used throughout the script
rng = np.random.default_rng(rng_seed)

# (Corrected and original) data from Hausser et al. (2019) are imported, as well as relative divergence data
yeast_df_name = 'rates_Hausser.csv'
yeast_df_path = os.path.join(path_data, yeast_df_name)
yeast_df = pd.read_csv(yeast_df_path)

folds_df_name = 'couples_divergence.csv'
folds_df_path = os.path.join(path_data, folds_df_name)
folds_df = pd.read_csv(folds_df_path)

# The appropriate subset of the data is selected
if anc_dist == 'all':
    yeast_data = yeast_df

elif anc_dist == 'dupli':
    yeast_data = yeast_df[(yeast_df['Duplication'] == 'SSD')|(yeast_df['Duplication'] == 'WGD')].reset_index(drop=True)

elif anc_dist == 'single':
    yeast_data = yeast_df[yeast_df['Duplication'] == 'S'].reset_index(drop=True)

if dataset == 'Original':
    bp_param = 'bp'
    bp_fold_param = 'bp_fold_original'
    pEst_fold_param = 'pEst_fold_original'


elif dataset == 'Corrected':
    bp_param = 'bp_calc'
    bp_fold_param = 'bp_fold_calc'
    pEst_fold_param = 'pEst_fold_calc'

# Minimal and maximal values for transcription and translation rates are set according to the specified parameters
if bm_bounds:
    bm_lims = ARGS.bm_bounds[1:-1].split(', ')

    if bm_lims[0] == 'None':
        bm_min = 0

    else:
        bm_min = 10**float(bm_lims[0])

    if bm_lims[1] == 'None':
        bm_max = np.max(10 ** yeast_data['bm_calc'])

    else:
        bm_max = 10**float(bm_lims[1])

else:
    bm_min = 0
    bm_max = np.max(10**yeast_data['bm_calc'])

if bp_bounds:
    bp_lims = ARGS.bp_bounds[1:-1].split(', ')

    if bp_lims[0] == 'None':
        bp_min = 0

    else:
        bp_min = 10**float(bp_lims[0])

    if bp_lims[1] == 'None':
        bp_max = np.max(10**yeast_data[bp_param])

    else:
        bp_max = 10**float(bp_lims[1])

else:
    bp_min = 0
    bp_max = np.max(10**yeast_data[bp_param])

# F) Mutational effects distributions
bivariate = ARGS.bivariate

ratio = ARGS.mut_ratio
p_bm = ratio/(1 + ratio)

bm_params = ARGS.bm_params[1:-1].split(', ')
bp_params = ARGS.bp_params[1:-1].split(', ')
skew_param = ARGS.skew

# If a bivariate distribution of mutational effects is used, the mut_ratio is used to set the two standard deviations
# accordingly. A brute-force optimization is performed to find the standard deviations of transcriptional and
# translational effects that give the desired mutational bias while keeping the mean absolute protein abundance change
# at the desired level.

if bivariate:

    # Means of transcriptional and translational mutational effects
    mu_bm = float(bm_params[0])
    mu_bp = float(bp_params[0])

    # Standard deviations of transcriptional and mutational effects (obtained through the brute-force optimization)
    sd_bm = float(bm_params[1])
    sd_bp = float(bp_params[1])

    sd_bp_opt = evol_funct.optimize_diff((mu_bm, mu_bp), ratio, correlation, (sd_bm, sd_bp), skew_param)
    sd_bm_opt = sd_bp_opt * ratio

    sigma_bm = sd_bm_opt
    sigma_bp = sd_bp_opt

    # The optimized SDs are kept in a file, so that it's possible to validate the optimization later
    sd_data = [sd_bm, sd_bp, sd_bm_opt, sd_bp_opt, ratio, run_name]
    sd_final = pd.DataFrame([sd_data,], columns=['Specified_SD_Bm', 'Specified_SD_Bp', 'Optimized_SD_Bm',
                                                 'Optimized_SD_Bp', 'Mut_ratio', 'run'])
    sd_final.to_csv(os.path.join(path_folder, f'SD_opt{run_name}.csv'))

else:
    # Transcription rate mutational effects
    mu_bm = float(bm_params[0])
    sigma_bm = float(bm_params[1])

    # Translation rate mutational effects
    mu_bp = float(bp_params[0])
    sigma_bp = float(bp_params[1])

# Covariance matrix if a bivariate normal distribution is used
if bivariate:
    cov = correlation * math.sqrt(sigma_bm**2 * sigma_bp**2)
    cov_mat = np.array([[sigma_bm**2, cov], [cov, sigma_bp**2]])

# 2- Definition of functions
# Choice of fixation probability function

fix_prob = evol_funct.metropolis

# 3- Initialization of the simulation
# Ancestral singletons are sampled according to the specified sampling method

# Constants are set to be able to test whether the sampled fitness functions are steep enough for the immediate
# loss of a paralog to be deleterious following the duplication event

# We first need to set the variation of the optimum of cumulative protein abundance following the duplication
total_bm = ARGS.total_bm

opt_change = dupli_effect
# Interesting choices: between 1.87 and 2.15, to avoid negative fitness values immediately after duplication

# Then, the test of fitness function curvature can be set up
Ne_scaling = total_bm * (1 - ((3 * total_bm) / (4 * opt_change)))

# First, if all parameters are randomized:
# Before, distributions of lm and Q (which may be used outside of this loop) are generated:
if lm_val:
    lm = lm_val
    lm_tocalc = lm_val

else:
    lm_tocalc = 10 ** yeast_data['lm']
    lm = rng.choice(lm_tocalc, size=n_couples)

if q_val:
    Q_dist = q_val

else:
    Q_dist = (c_m * alpha_m * lm_tocalc) / (10 ** yeast_data[bp_param] / 10 ** yeast_data['bm_calc'])

# Only 21 genes have a higher noise sensitivity that Qmax (6.8588e-06) predicted by Hausser et al.
# They are dropped from the Q_dist distribution prior to sampling.
Q_dist = Q_dist[Q_dist <= 6.8588e-06]
Q_sim = rng.choice(Q_dist, size=n_couples)

# Then, the sampling using all-randomized parameters can be performed
if samp_meth == 'random':

    # Distribution of optimal protein abundances
    pOpt_dist = (10 ** yeast_data['bm_calc'] * 10 ** yeast_data[bp_param]) / (alpha_m * alpha_p)
    pOpt_sim = rng.choice(pOpt_dist, size=n_couples)

    # The next steps select only genes for which the fitness function is steep enough for the immediate loss of a
    # paralog to be deleterious. This is ignored in the case that a Q value has been specified manually (presumably to
    # perform a sensitivity analysis according to this parameter).

    if not q_val:
        # Only ancestral genes for which the loss of a duplicate will be deleterious immediately after the duplication
        # are kept.
        # Sampling is done again for the genes for which the fitness function is not steep enough
        curve_test = np.where(Q_sim * pOpt_sim <= 1 / (Ne_scaling * pop_size), True, False)

        while np.any(curve_test):
            nums = np.where(curve_test)

            for gene in nums[0]:
                Q_sim[gene] = rng.choice(Q_dist, size=1)
                pOpt_sim[gene] = rng.choice(pOpt_dist, size=1)
            curve_test = np.where(Q_sim * pOpt_sim <= 1 / (Ne_scaling * pop_size), True, False)

    # The resulting singletons then need to be assigned optimal Bm and Bp rates according to the cost-precision
    # trade-off (Hausser et al., 2019). Adequate ancestral states also need to be generated for the four "parallel"
    # simulations

    anc_mixed = pd.DataFrame(columns=['pOpt', 'Q', 'lm', 'Bm', 'Bp'])
    anc_mixed['pOpt'] = pOpt_sim
    anc_mixed['Q'] = Q_sim
    anc_mixed['lm'] = lm

    anc_ADS = anc_mixed.copy()

    # The optimization of singletons according to the cost-precision trade-off is first performed. The resulting
    # ancestral states will be used for both the "mixed" and "no-cost" simulations
    print('Optimizing the ancestral singletons.')

    for row in range(anc_mixed.shape[0]):
        Q = anc_mixed.at[row, 'Q']
        pOpt = anc_mixed.at[row, 'pOpt']
        length = anc_mixed.at[row, 'lm']

        bm_est = math.sqrt(pOpt * alpha_p * alpha_m * (Q / (c_m * length * alpha_m)))  # Bm according to Hausser et al.
        bm_bounds = ((bm_est - (0.5 * bm_est)), (bm_est + (0.5 * bm_est)))

        bp_est = (pOpt * alpha_m * alpha_p) / bm_est
        bp_bounds = ((bp_est - (0.5 * bp_est)), (bp_est + (0.5 * bp_est)))

        opt_single = optimize.differential_evolution(evol_funct.ancestral_opt, (bm_bounds, bp_bounds), tol=1e-10,
                                                     args=(pOpt, Q, length, c_m, alpha_m, alpha_p, cv_0), seed=rng_seed)

        bm_optimized = opt_single.x[0]
        bp_optimized = opt_single.x[1]

        anc_mixed.at[row, 'Bm'] = bm_optimized
        anc_mixed.at[row, 'Bp'] = bp_optimized

        if row % 100 == 0:
            print(f'Done with gene {row}')

    print('Optimization of the ancestral states complete.')

    # Some of the resulting optimal singleton genes might be outside the previously defined boundaries for Bm and Bp.
    # Because the optimization process is a bit tedious, any invalid gene is simply filtered out and the remaining
    # singletons are sampled with replacement to obtain the desired number of gene pairs.
    anc_mixed = anc_mixed[(anc_mixed['Bm'] > bm_min) & (anc_mixed['Bp'] > bp_min)].reset_index(drop=True)
    anc_mixed = anc_mixed[(anc_mixed['Bm'] < bm_max) & (anc_mixed['Bp'] < bp_max)].reset_index(drop=True)

    n_missing = n_couples - anc_mixed.shape[0]
    n_first = n_missing

    if n_missing != 0:
        to_add = rng.choice(anc_mixed, size=n_missing)
        to_add = pd.DataFrame(to_add, columns=anc_mixed.columns)

        anc_mixed = pd.concat([anc_mixed, to_add]).reset_index(drop=True)

        print(f'{n_first} invalid singletons were replaced (cost-precision).')

    # Ancestral states for the ADS-only simulation are next defined. In that case, the optimal Bm and Bp according to
    # Hausser et al. (2019) are simply used as the ancestral optima
    anc_ADS['Bm'] = np.sqrt(anc_ADS['pOpt'] * alpha_p * alpha_m * (anc_ADS['Q'] / (c_m * anc_ADS['lm'] * alpha_m)))
    anc_ADS['Bp'] = (anc_ADS['pOpt'] * alpha_m * alpha_p) / anc_ADS['Bm']

    # In case some genes are outside the specified boundaries for Bm and Bp, the same approach is used as previously:
    anc_ADS = anc_ADS[(anc_ADS['Bm'] > bm_min) & (anc_ADS['Bp'] > bp_min)].reset_index(drop=True)
    anc_ADS = anc_ADS[(anc_ADS['Bm'] < bm_max) & (anc_ADS['Bp'] < bp_max)].reset_index(drop=True)

    n_missing = n_couples - anc_ADS.shape[0]
    n_first = n_missing

    if n_missing != 0:
        to_add = rng.choice(anc_ADS, size=n_missing)
        to_add = pd.DataFrame(to_add, columns=anc_ADS.columns)

        anc_ADS = pd.concat([anc_ADS, to_add]).reset_index(drop=True)

        print(f'{n_first} invalid singletons were replaced (ADS-only).')

    # Ancestral states are generated for a fourth "minimal" simulation. It is also strictly ADS, but this time noise
    # sensitivities (Q) are randomized to avoid any influence of the cost-precision constraints on the distribution of
    # ancestral genes in the expression space.
    anc_min = anc_ADS.copy()
    anc_min['Q'] = rng.choice(Q_dist, size=n_couples)

    # As before, a sufficient curvature of the fitness function is ensured
    curve_test = np.where(anc_min['Q'] * anc_min['pOpt'] <= 1 / (Ne_scaling * pop_size), True, False)

    while np.any(curve_test):
        nums = np.where(curve_test)

        for gene in nums[0]:
            anc_min.at[gene, 'Q'] = rng.choice(Q_dist, size=1)
        curve_test = np.where(anc_min['Q'] * anc_min['pOpt'] <= 1 / (Ne_scaling * pop_size), True, False)

    # As this dataframe was generated from anc_ADS, any invalid Bm-Bp combinations have already been filtered out
    # and the generation of the ancestral states is thus complete.

# Second, if actual genes are sampled
elif samp_meth == 'genes':
    # A dataframe containing only genes with a steep enough fitness function (to ensure that the immediate loss of a
    # newly created paralog would be deleterious) is first created
    yeast_tosamp = yeast_data[['bm_calc', bp_param, 'lm']].copy()
    yeast_tosamp['pOpt'] = (10**yeast_tosamp['bm_calc'] * 10**yeast_tosamp[bp_param]) / (alpha_m * alpha_p)
    yeast_tosamp['Q'] = (c_m * alpha_m * 10**yeast_tosamp['lm']) / (10**yeast_tosamp[bp_param] / 10**yeast_tosamp['bm_calc'])
    yeast_tosamp['curve_test'] = yeast_tosamp['Q'] * yeast_tosamp['pOpt']
    yeast_tosamp = yeast_tosamp[yeast_tosamp['curve_test'] >= 1 / (Ne_scaling * pop_size)].reset_index(drop=True)

    yeast_tosamp['bm_calc'] = 10**yeast_tosamp['bm_calc']
    yeast_tosamp[bp_param] = 10**yeast_tosamp[bp_param]
    yeast_tosamp['lm'] = 10**yeast_tosamp['lm']

    # Then, a set of genes can be sampled
    anc_mixed = pd.DataFrame(columns=['pOpt', 'Q', 'lm', 'Bm', 'Bp'])

    gen_sample = rng.choice(yeast_tosamp[['bm_calc', bp_param, 'lm', 'pOpt', 'Q']], size=n_couples)
    gen_sample = pd.DataFrame(gen_sample, columns=['bm_calc', bp_param, 'lm', 'pOpt', 'Q'])

    # Any gene outside the specified boundaries for Bm and Bp is replaced before any further manipulations
    gen_sample = gen_sample[(gen_sample['bm_calc'] > bm_min) & (gen_sample[bp_param] > bp_min)]
    gen_sample = gen_sample[(gen_sample['bm_calc'] < bm_max) & (gen_sample[bp_param] < bp_max)].reset_index(drop=True)

    n_missing = n_couples - gen_sample.shape[0]

    while n_missing != 0:
        to_add = rng.choice(yeast_tosamp[['bm_calc', bp_param, 'lm', 'pOpt', 'Q']], size=n_missing)
        to_add = pd.DataFrame(to_add, columns=['bm_calc', bp_param, 'lm', 'pOpt', 'Q'])
        gen_sample = pd.concat([gen_sample, to_add]).reset_index(drop=True)

        gen_sample = gen_sample[(gen_sample['bm_calc'] > bm_min) & (gen_sample[bp_param] > bp_min)]
        gen_sample = gen_sample[(gen_sample['bm_calc'] < bm_max) & (gen_sample[bp_param] < bp_max)].reset_index(drop=True)
        n_missing = n_couples - gen_sample.shape[0]

    # The anc_mixed dataframe is filled
    anc_mixed['Bm'] = gen_sample['bm_calc']
    anc_mixed['Bp'] = gen_sample[bp_param]
    anc_mixed['lm'] = gen_sample['lm']
    anc_mixed['pOpt'] = gen_sample['pOpt']
    anc_mixed['Q'] = gen_sample['Q']

    # For the mixed and no_cost simulations, these newly sampled genes are optimized. This should only slightly change
    # their Bm and Bp rates. This is still followed by a filtering step to unsure that all ancestral genes are within
    # the defined boundaries for Bm and Bp.
    anc_ADS = anc_mixed.copy()

    print('Optimizing the ancestral singletons.')

    for row in range(anc_mixed.shape[0]):
        Q = anc_mixed.at[row, 'Q']
        pOpt = anc_mixed.at[row, 'pOpt']
        length = anc_mixed.at[row, 'lm']

        bm_est = math.sqrt(pOpt * alpha_p * alpha_m * (Q / (c_m * length * alpha_m)))  # Bm according to Hausser et al.
        bm_bounds = ((bm_est - (0.5 * bm_est)), (bm_est + (0.5 * bm_est)))

        bp_est = (pOpt * alpha_m * alpha_p) / bm_est
        bp_bounds = ((bp_est - (0.5 * bp_est)), (bp_est + (0.5 * bp_est)))

        opt_single = optimize.differential_evolution(evol_funct.ancestral_opt, (bm_bounds, bp_bounds), tol=1e-10,
                                                     args=(pOpt, Q, length, c_m, alpha_m, alpha_p, cv_0), seed=rng_seed)

        bm_optimized = opt_single.x[0]
        bp_optimized = opt_single.x[1]

        anc_mixed.at[row, 'Bm'] = bm_optimized
        anc_mixed.at[row, 'Bp'] = bp_optimized

        if row % 100 == 0:
            print(f'Done with gene {row}')

    print('Optimization of the ancestral states complete.')

    # Filtering of any genes made invalid by this optimization step
    anc_mixed = anc_mixed[(anc_mixed['Bm'] > bm_min) & (anc_mixed['Bp'] > bp_min)].reset_index(drop=True)
    anc_mixed = anc_mixed[(anc_mixed['Bm'] < bm_max) & (anc_mixed['Bp'] < bp_max)].reset_index(drop=True)

    n_missing = n_couples - anc_mixed.shape[0]
    n_first = n_missing

    if n_missing != 0:
        to_add = rng.choice(anc_mixed, size=n_missing)
        to_add = pd.DataFrame(to_add, columns=anc_mixed.columns)

        anc_mixed = pd.concat([anc_mixed, to_add]).reset_index(drop=True)

        print(f'{n_first} invalid singletons were replaced (cost-precision).')

    # Ancestral states for the ADS-only simulation are next defined, as previously.
    anc_ADS['Bm'] = np.sqrt(anc_ADS['pOpt'] * alpha_p * alpha_m * (anc_ADS['Q'] / (c_m * anc_ADS['lm'] * alpha_m)))
    anc_ADS['Bp'] = (anc_ADS['pOpt'] * alpha_m * alpha_p) / anc_ADS['Bm']

    # In case some genes are outside the specified boundaries for Bm and Bp, the same filtering approach is used:
    anc_ADS = anc_ADS[(anc_ADS['Bm'] > bm_min) & (anc_ADS['Bp'] > bp_min)].reset_index(drop=True)
    anc_ADS = anc_ADS[(anc_ADS['Bm'] < bm_max) & (anc_ADS['Bp'] < bp_max)].reset_index(drop=True)

    n_missing = n_couples - anc_ADS.shape[0]
    n_first = n_missing

    if n_missing != 0:
        to_add = rng.choice(anc_ADS, size=n_missing)
        to_add = pd.DataFrame(to_add, columns=anc_ADS.columns)

        anc_ADS = pd.concat([anc_ADS, to_add]).reset_index(drop=True)

        print(f'{n_first} invalid singletons were replaced (ADS-only).')

    # Finally, the ancestral states are generated for the fourth "minimal" simulation.
    anc_min = anc_ADS.copy()
    anc_min['Q'] = rng.choice(Q_dist, size=n_couples)

    # As before, a sufficient curvature of the fitness function is ensured
    curve_test = np.where(anc_min['Q'] * anc_min['pOpt'] <= 1 / (Ne_scaling * pop_size), True, False)

    while np.any(curve_test):
        nums = np.where(curve_test)

        for gene in nums[0]:
            anc_min.at[gene, 'Q'] = rng.choice(Q_dist, size=1)
        curve_test = np.where(anc_min['Q'] * anc_min['pOpt'] <= 1 / (Ne_scaling * pop_size), True, False)

    # As this dataframe was generated from anc_ADS, any invalid Bm-Bp combinations have already been filtered out
    # and the generation of the ancestral states is thus complete.

# The initial post-duplication states are then generated from the ancestral singletons obtained above
pOpt_dupli_mixed = anc_mixed['pOpt'] * opt_change  # Also for the no-cost simulation
pOpt_dupli_ADS = anc_ADS['pOpt'] * opt_change
pOpt_dupli_min = anc_min['pOpt'] * opt_change

# Post-duplication transcription rates
post_bm_mixed = (total_bm * anc_mixed['Bm']) / 2  # Also for the no-cost simulation
post_bm_ADS = (total_bm * anc_ADS['Bm']) / 2
post_bm_min = (total_bm * anc_min['Bm'] / 2)

rates_init = np.zeros((n_couples, 8))
rates_init[:, 0] = range(n_couples)
rates_init[:, 1] = anc_mixed['Q']  # Q for each gene pair
rates_init[:, 2] = pOpt_dupli_mixed  # post-duplication pOpt for each gene pair
rates_init[:, 3] = anc_mixed['lm']  # Length of the mRNAs (important in case it varies between paralog pairs)
rates_init[:, 4] = post_bm_mixed
rates_init[:, 6] = post_bm_mixed
rates_init[:, 5] = anc_mixed['Bp']
rates_init[:, 7] = anc_mixed['Bp']
rates_mixed = rates_init.copy()
rates_nocost = rates_init.copy()  # To also run the simulation while neglecting the cost of transcription

# For the ADS-only simulation
start_ADS = np.zeros((n_couples, 8))
start_ADS[:, 0] = range(n_couples)
start_ADS[:, 1] = anc_ADS['Q']  # Q for each gene pair
start_ADS[:, 2] = pOpt_dupli_ADS  # post-duplication pOpt for each gene pair
start_ADS[:, 3] = anc_ADS['lm']  # Length of the mRNAs (important in case it varies between paralog pairs)
start_ADS[:, 4] = post_bm_ADS
start_ADS[:, 6] = post_bm_ADS
start_ADS[:, 5] = anc_ADS['Bp']
start_ADS[:, 7] = anc_ADS['Bp']
rates_ADS = start_ADS.copy()

# For the minimal simulation
start_min = np.zeros((n_couples, 8))
start_min[:, 0] = range(n_couples)
start_min[:, 1] = anc_min['Q']  # Q for each gene pair
start_min[:, 2] = pOpt_dupli_min  # post-duplication pOpt for each gene pair
start_min[:, 3] = anc_ADS['lm']  # Length of the mRNAs (important in case it varies between paralog pairs)
start_min[:, 4] = post_bm_min
start_min[:, 6] = post_bm_min
start_min[:, 5] = anc_min['Bp']
start_min[:, 7] = anc_min['Bp']
rates_min = start_min.copy()

# Then, the data structures used to save expression rates at the beginning and the end of the simulation are made.

# One dataframe per selection regime. Each row is one duplicate couple at a specific mutation round, and columns are
# as follows: Round, Couple, Q, pOpt, lm, Bm1, Bm2, Bp1, Bp2, Prot1, Prot2, cv1, cv2, Exp_cost

data_model = pd.DataFrame(columns=['Round', 'Couple', 'Q', 'pOpt', 'lm', 'Bm1', 'Bp1', 'Bm2', 'Bp2', 'Prot1', 'Prot2',
                                   'cv1', 'cv2', 'Exp_cost', 'Fitness'])
data_model['Couple'] = range(n_couples)
data_model['Round'] = 0

init_noise = data_model.copy()
init_noise.iloc[:, 2:9] = rates_init[:, 1:8]
init_noise.iloc[:, 9] = ((init_noise['Bm1'] * init_noise['Bp1']) / (alpha_m * alpha_p))
init_noise.iloc[:, 10] = ((init_noise['Bm2'] * init_noise['Bp2']) / (alpha_m * alpha_p))
init_noise.iloc[:, 11] = np.sqrt((init_noise['Prot1']**2)*(1/init_noise['Prot1']) + (alpha_p/init_noise['Bm1']) + cv_0**2)
init_noise.iloc[:, 12] = np.sqrt((init_noise['Prot2']**2)*(1/init_noise['Prot2']) + (alpha_p/init_noise['Bm2']) + cv_0**2)
init_noise['Exp_cost'] = (init_noise['Bm1'] + init_noise['Bm2']) * lm * c_m

init_mixed = init_noise.copy()
init_nocost = init_noise.copy()

# Addition of fitness values
init_mixed['Fitness'] = evol_funct.fit_global_dupli(rates_mixed[:, 4], rates_mixed[:, 6], rates_mixed[:, 5],
                                                    rates_mixed[:, 7], rates_mixed[:, 2], rates_mixed[:, 1], alpha_m,
                                                    alpha_p, cv_0, rates_mixed[:, 3], c_m)
init_nocost['Fitness'] = evol_funct.fit_noise_dupli(rates_nocost[:, 4], rates_nocost[:, 6], rates_nocost[:, 5],
                                                    rates_nocost[:, 7], rates_nocost[:, 2], rates_nocost[:, 1], alpha_m,
                                                    alpha_p, cv_0)

init_ADS = data_model.copy()
init_ADS.iloc[:, 2:9] = rates_ADS[:, 1:8]
init_ADS.iloc[:, 9] = ((init_ADS['Bm1'] * init_ADS['Bp1']) / (alpha_m * alpha_p))
init_ADS.iloc[:, 10] = ((init_ADS['Bm2'] * init_ADS['Bp2']) / (alpha_m * alpha_p))
init_ADS.iloc[:, 11] = np.sqrt((init_ADS['Prot1']**2)*(1/init_ADS['Prot1']) + (alpha_p/init_ADS['Bm1']) + cv_0**2)
init_ADS.iloc[:, 12] = np.sqrt((init_ADS['Prot2']**2)*(1/init_ADS['Prot2']) + (alpha_p/init_ADS['Bm2']) + cv_0**2)
init_ADS['Exp_cost'] = (init_ADS['Bm1'] + init_ADS['Bm2']) * lm * c_m

init_min = data_model.copy()
init_min.iloc[:, 2:9] = rates_min[:, 1:8]
init_min.iloc[:, 9] = ((init_min['Bm1'] * init_min['Bp1']) / (alpha_m * alpha_p))
init_min.iloc[:, 10] = ((init_min['Bm2'] * init_min['Bp2']) / (alpha_m * alpha_p))
init_min.iloc[:, 11] = np.sqrt((init_min['Prot1']**2)*(1/init_min['Prot1']) + (alpha_p/init_min['Bm1']) + cv_0**2)
init_min.iloc[:, 12] = np.sqrt((init_min['Prot2']**2)*(1/init_min['Prot2']) + (alpha_p/init_min['Bm2']) + cv_0**2)
init_min['Exp_cost'] = (init_min['Bm1'] + init_min['Bm2']) * lm * c_m

# Addition of fitness values
init_ADS['Fitness'] = evol_funct.fit_parabola(rates_ADS[:, 4], rates_ADS[:, 6], rates_ADS[:, 5], rates_ADS[:, 7],
                                              rates_ADS[:, 2], rates_ADS[:, 1], alpha_m, alpha_p)
init_min['Fitness'] = evol_funct.fit_parabola(rates_min[:, 4], rates_min[:, 6], rates_min[:, 5], rates_min[:, 7],
                                              rates_min[:, 2], rates_min[:, 1], alpha_m, alpha_p)

# Combining all four initial dataframes and saving them in csv format
init_mixed['Model'] = 'Mixed'
init_nocost['Model'] = 'No Cost'
init_ADS['Model'] = 'ADS'
init_min['Model'] = 'Minimal'
data_all = pd.concat([init_mixed, init_nocost, init_ADS, init_min]).reset_index(drop=True)
#init_all.to_csv(os.path.join(path_folder, f'data_all_{run_name}.csv'), index=False)

# Data structures are initialized to keep track of the number of fixed mutations through time
muta_model = pd.DataFrame(columns=['Round', 'Couple', 'P1 Mutations', 'P2 Mutations'])
muta_model['Couple'] = range(n_couples)
muta_model['P1 Mutations'] = 0
muta_model['P2 Mutations'] = 0

muta_mixed = muta_model.copy()
muta_mixed['Round'] = 0
muta_nocost = muta_mixed.copy()
muta_ADS = muta_mixed.copy()
muta_min = muta_mixed.copy()

muta_mixed_current = muta_model.copy()
muta_nocost_current = muta_model.copy()
muta_ADS_current = muta_model.copy()
muta_min_current = muta_model.copy()

# Another dataframe is generated (and saved as csv). It will be used to store the p-values of a Mood's test assessing
# whether the median of the distribution of protein abundance log2 fold-changes associated with each selection regimes
# is significantly different from the "true" distribution for yeast duplicates. This is followed throughout the
# mutation rounds, and can be used as a end condition for the simulation.

# The "true" distribution of log2 protein abundance fold-changes calculated from Hausser's data (after correction)
# is imported
if dupli_type == 'all':
    fold_med = folds_df.copy()

elif dupli_type == 'WGD':
    fold_med = folds_df[folds_df['Duplication'] == 'WGD'].copy()

elif dupli_type == 'SSD':
    fold_med = folds_df[folds_df['Duplication'] == 'SSD'].copy()

# The dataframe of p-values is initialized
Mood_pvals_model = pd.DataFrame(columns=['Round', 'Mixed', 'No Cost', 'ADS Only', 'Minimal'])
Mood_pvals_model.at[0, 'Round'] = 0
Mood_pvals = Mood_pvals_model.copy()

fold_init = evol_funct.fold_change('Prot1', 'Prot2', init_noise)
fold_init_ADS = evol_funct.fold_change('Prot1', 'Prot2', init_ADS)

pval_init = stats.median_test(fold_med[pEst_fold_param], fold_init)[1]
pval_init_ADS = stats.median_test(fold_med[pEst_fold_param], fold_init_ADS)[1]

Mood_pvals.at[0, 'Mixed'] = pval_init
Mood_pvals.at[0, 'No Cost'] = pval_init
Mood_pvals.at[0, 'ADS Only'] = pval_init_ADS
Mood_pvals.at[0, 'Minimal'] = pval_init_ADS

# If gene loss is allowed, a data structure is initialized to keep track or the number of duplicate couples through time
if dupli_loss:
    loss_df = pd.DataFrame(columns=['Round', 'Mixed', 'No Cost', 'ADS Only', 'Minimal'])
    loss_df.at[0, 'Round'] = 0

    loss_df.iloc[0, 1:4] = n_couples
    loss_df.to_csv(os.path.join(path_folder, f'Gene_loss_{run_name}.csv'), index=False)

# Before the simulation is started, one figure is generated to validate the optimality of the ancestral singleton states
# The figure is constructed for 50 couples selected at random, and all figures are saved in the same multi-pages pdf
if make_figs:
    ran_couples = rng.choice(range(n_couples), size=50)
    anc_fig = PdfPages(os.path.join(path_folder, f'Ancestral_states_{run_name}.pdf'))

    for couple in ran_couples:
        bm_mixed = anc_mixed.at[couple, 'Bm']
        bp_mixed = anc_mixed.at[couple, 'Bp']
        bm_ADS = anc_ADS.at[couple, 'Bm']
        bp_ADS = anc_ADS.at[couple, 'Bp']
        pOpt_mixed = anc_mixed.at[couple, 'pOpt']
        pOpt_ADS = anc_ADS.at[couple, 'pOpt']
        Q_mixed = anc_mixed.at[couple, 'Q']
        Q_ADS = anc_ADS.at[couple, 'Q']
        lm = anc_mixed.at[couple, 'lm']

        panel_couple = evol_funct.panel_ancestral(bm_mixed, bp_mixed, bm_ADS, bp_ADS, pOpt_mixed, pOpt_ADS, Q_mixed,
                                                  Q_ADS, alpha_m, alpha_p, cv_0, lm, c_m)

        panel_couple.savefig(anc_fig, format='pdf')
        plt.close()

    anc_fig.close()

# The simulations can then be performed:
stop_sim = False   # To end the simulation at some point

done_mixed = False  # To check whether the end condition has been met in any of the four simulations
done_nocost = False
done_ADS = False
done_min = False

last_mixed = False  # To check if the last mutation round of a given simulation has just been made
last_nocost = False
last_ADS = False
last_min = False

step = 1

print(f'Starting run {run_name}.')

# Reference time for the mutation-selection rounds only
t0_evol = time.process_time()

while not stop_sim:
    # The ancestral (pre-mutation) fitness is first computed for each selection regime
    fit_mixed = evol_funct.fit_global_dupli(rates_mixed[:, 4], rates_mixed[:, 6], rates_mixed[:, 5],
                                            rates_mixed[:, 7], rates_mixed[:, 2], rates_mixed[:, 1], alpha_m, alpha_p,
                                            cv_0, rates_mixed[:, 3], c_m)
    fit_nocost = evol_funct.fit_noise_dupli(rates_nocost[:, 4], rates_nocost[:, 6], rates_nocost[:, 5],
                                            rates_nocost[:, 7], rates_nocost[:, 2], rates_nocost[:, 1],
                                            alpha_m, alpha_p, cv_0)
    fit_ADS = evol_funct.fit_parabola(rates_ADS[:, 4], rates_ADS[:, 6], rates_ADS[:, 5], rates_ADS[:, 7],
                                      rates_ADS[:, 2], rates_ADS[:, 1], alpha_m, alpha_p)
    fit_min = evol_funct.fit_parabola(rates_min[:, 4], rates_min[:, 6], rates_min[:, 5], rates_min[:, 7],
                                      rates_min[:, 2], rates_min[:, 1], alpha_m, alpha_p)

    # Mutations are generated, either using a bivariate distribution or randomly selecting independent
    # transcriptional and translational effects
    if not bivariate:

        # For each couple, one of the two rates is chosen to be mutated, according to the level of mutational bias
        bm_decision = rng.random(size=n_couples)
        bm_choice = np.where(bm_decision <= p_bm, 1, 0)
        bp_choice = np.abs(bm_choice - 1)

        # For the case where mutations independently affect either transcription or translation, two approaches are
        # possible. The first considers that any asymmetry of mutational effects is unidirectional
        if not bidirectional:
            # An array of bm and bp mutations (identical for each selection regime) is generated
            bm_muta = stats.skewnorm.rvs(skew_param, loc=mu_bm, scale=sigma_bm, size=n_couples, random_state=rng)
            bp_muta = stats.skewnorm.rvs(skew_param, loc=mu_bp, scale=sigma_bp, size=n_couples, random_state=rng)

            # Choices generated before are applied
            bm_muta = bm_muta * bm_choice
            bp_muta = bp_muta * bp_choice

            # Matrix of mutational effects
            mutations = np.zeros((n_couples, 4))
            mutations[:, 0] = bm_muta
            mutations[:, 1] = bp_muta
            mutations[:, 2] = bm_muta
            mutations[:, 3] = bp_muta

        # The second considers that high transcription/translation rates are biased towards a decrease, while low
        # rates are biased towards an increase through mutations
        elif bidirectional:
            skew_abs = abs(skew_param)
            bm_to_low = stats.skewnorm.rvs(-skew_abs, loc=mu_bm, scale=sigma_bm, size=n_couples,
                                           random_state=rng) * bm_choice
            bm_to_high = stats.skewnorm.rvs(skew_abs, loc=mu_bm, scale=sigma_bm, size=n_couples,
                                            random_state=rng) * bm_choice
            bm_mid = stats.skewnorm.rvs(0, loc=mu_bm, scale=sigma_bm, size=n_couples, random_state=rng) * bm_choice

            bp_to_low = stats.skewnorm.rvs(-skew_abs, loc=mu_bp, scale=sigma_bp, size=n_couples,
                                           random_state=rng) * bp_choice
            bp_to_high = stats.skewnorm.rvs(skew_abs, loc=mu_bp, scale=sigma_bp, size=n_couples,
                                            random_state=rng) * bp_choice
            bp_mid = stats.skewnorm.rvs(0, loc=mu_bm, scale=sigma_bp, size=n_couples, random_state=rng) * bp_choice

            # The thresholds are set
            low_lim_bm = np.quantile(yeast_data['bm_calc'], ext_quant)
            high_lim_bm = np.quantile(yeast_data['bm_calc'], (1 - ext_quant))
            low_lim_bp = np.quantile(yeast_data[bp_param], ext_quant)
            high_lim_bp = np.quantile(yeast_data[bp_param], (1 - ext_quant))

    elif bivariate:
        # An array of mutations, each with transcriptional and translational effects, is generated
        corr_muta = evol_funct.multivariate_skewnorm((skew_param, skew_param), (mu_bm, mu_bp),
                                                     cov=cov_mat).rvs_fast(size=n_couples, random_state=rng)
        # Matrix of mutational effects
        mutations = np.zeros((n_couples, 4))
        mutations[:, 0:2] = corr_muta
        mutations[:, 2:4] = corr_muta

    # Then, copy P1 or P2 is randomly chosen to receive the mutation in each couple
    # The same P1_choice array will be used later to generate the final array of mutations if bidirectional asymmetry
    # has been enabled
    P1_decision = rng.random(size=n_couples)
    P1_choice = np.where(P1_decision <= 0.5, 1, 0)
    P2_choice = np.abs(P1_choice - 1)

    # If the potential asymmetry of the mutational effects distributions has not been set to be bidirectional, the final
    # arrays of mutations can be generated. Otherwise, it is done later, separately for each of the four simulations.
    if not bidirectional:
        mutations[:, 0:2] = mutations[:, 0:2] * P1_choice[:, np.newaxis]
        mutations[:, 2:4] = mutations[:, 2:4] * P2_choice[:, np.newaxis]

    # Mutational effects are applied, fixation probabilities are computed, and mutations reach fixation or are lost.
    # This is done separately for each selection regime, depending on whether the end condition has already been
    # reached. The same set of random numbers are used in the four cases to decide if mutations fix (according to pfix).

    fix_decision = rng.random(size=n_couples)

    if dupli_loss:
        # Preparing the dataframe to save the number of intact couples through time
        loss_data = pd.DataFrame(columns=['Round', 'Mixed', 'No Cost', 'ADS', 'Minimal'])

        loss_data.at[0, 'Round'] = step

    # First, for the mixed model
    if not done_mixed:

        if dupli_loss:
            # Gene loss is performed when neutral
            loss_mixed = evol_funct.gene_loss(rates_mixed, fit_mixed, evol_funct.fit_global_dupli, pop_size,
                                              args_fit=(alpha_m, alpha_p, cv_0, rates_mixed[:, 3], c_m))

            # The rates dataframe is updated to take into account the newly performed gene loss
            rates_mixed = loss_mixed[0].copy()

            # The remaining number of intact couples is saved
            n_mixed = loss_mixed[1]
            loss_data.at[0, 'Mixed'] = n_mixed

        # If bidirectional asymmetry of mutational effects has been used, the final array of mutations for this
        # simulation is generated
        if bidirectional:
            if bivariate:
                # Because this has not been implemented for a bivariate distribution of mutational effects
                raise Exception

            mutations_mixed = evol_funct.muta_bidirectional(rates_mixed, n_couples, (low_lim_bm, high_lim_bm),
                                                            (low_lim_bp, high_lim_bp), bm_to_low, bm_mid, bm_to_high,
                                                            bp_to_low, bp_mid, bp_to_high, P1_choice)

        # Otherwise, the array of mutations generated before is used as is
        else:
            mutations_mixed = mutations

        # The ancestral rates are copied ahead of the calculation of new mutations
        new_mixed = rates_mixed.copy()

        # Applying mutational effects. The ancestral rates are used to compute the absolute effects of mutations if
        # the epistasy has been set to additive
        if epistasy == 'mult':
            new_mixed[:, 4:8] = new_mixed[:, 4:8] + (new_mixed[:, 4:8] * mutations_mixed)

        elif epistasy == 'add':
            new_mixed[:, 4:8] = new_mixed[:, 4:8] + (rates_init[:, 4:8] * mutations_mixed)

        # Computing mutant fitness
        fit_new_mixed = evol_funct.fit_global_dupli(new_mixed[:, 4], new_mixed[:, 6], new_mixed[:, 5],
                                                    new_mixed[:, 7], new_mixed[:, 2], new_mixed[:, 1], alpha_m, alpha_p,
                                                    cv_0, new_mixed[:, 3], c_m)

        # Ancestral and mutant fitness are normalized between 0 and 1
        fit_mixed = fit_mixed / 0.42
        fit_new_mixed = fit_new_mixed / 0.42

        # Calculating fixation probability from the fitness difference (ancestral - mutant)
        #prob_mixed = np.where(fit_new_mixed > 0, fix_prob(fit_mixed, fit_new_mixed, pop_size), 0)
        prob_mixed = fix_prob(fit_mixed, fit_new_mixed, pop_size)

        # Decisions on mutations
        mixed_decision = np.where(fix_decision < prob_mixed, 1, 0)

        # If the mutation takes Bm or Bp below or above the boundaries, it is canceled
        low_bm1 = np.where(new_mixed[:, 4] < bm_min, 0, 1)
        low_bm2 = np.where(new_mixed[:, 6] < bm_min, 0, 1)
        low_bp1 = np.where(new_mixed[:, 5] < bp_min, 0, 1)
        low_bp2 = np.where(new_mixed[:, 7] < bp_min, 0, 1)
        mixed_decision = mixed_decision * (low_bm1 * low_bm2 * low_bp1 * low_bp2)

        high_bm1 = np.where(new_mixed[:, 4] > bm_max, 0, 1)
        high_bm2 = np.where(new_mixed[:, 6] > bm_max, 0, 1)
        high_bp1 = np.where(new_mixed[:, 5] > bp_max, 0, 1)
        high_bp2 = np.where(new_mixed[:, 7] > bp_max, 0, 1)
        mixed_decision = mixed_decision * (high_bm1 * high_bm2 * high_bp1 * high_bp2)

        # Identifying all couples which are not mutated in the current round
        kept_mixed = np.abs(mixed_decision - 1)

        # Counting new mutations in copies P1 and P2, and adding them to the previous totals
        mut_P1_mixed = P1_choice * mixed_decision
        mut_P2_mixed = P2_choice * mixed_decision

        muta_mixed_current['P1 Mutations'] += mut_P1_mixed
        muta_mixed_current['P2 Mutations'] += mut_P2_mixed

        # Creating the array where only the mutants are kept
        mutants_mixed = new_mixed.copy()
        mutants_mixed[:, 4:8] = mutants_mixed[:, 4:8] * mixed_decision[:, np.newaxis]

        # Creating the array where only the non-mutated couples are kept
        nomut_mixed = rates_mixed.copy()
        nomut_mixed[:, 4:8] = nomut_mixed[:, 4:8] * kept_mixed[:, np.newaxis]

        # Combining them in a final array
        final_mixed = mutants_mixed.copy()
        final_mixed[:, 4:8] = final_mixed[:, 4:8] + nomut_mixed[:, 4:8]

    else:
        final_mixed = rates_mixed.copy()

        if dupli_loss:
            loss_data.at[0, 'Mixed'] = n_mixed

    # Second, for the no-cost model
    if not done_nocost:

        if dupli_loss:
            # Gene loss is performed when neutral
            loss_nocost = evol_funct.gene_loss(rates_nocost, fit_nocost, evol_funct.fit_noise_dupli, pop_size,
                                               (alpha_m, alpha_p, cv_0))

            # The rates dataframe is updated to take into account the newly performed gene loss
            rates_nocost = loss_nocost[0].copy()

            # The remaining number of intact couples is saved
            n_nocost = loss_nocost[1]
            loss_data.at[0, 'No Cost'] = n_nocost

        # If bidirectional asymmetry of mutational effects has been used, the final array of mutations for this
        # simulation is generated
        if bidirectional:
            if bivariate:
                # Because this has not been implemented for a bivariate distribution of mutational effects
                raise Exception

            mutations_nocost = evol_funct.muta_bidirectional(rates_nocost, n_couples, (low_lim_bm, high_lim_bm),
                                                             (low_lim_bp, high_lim_bp), bm_to_low, bm_mid, bm_to_high,
                                                             bp_to_low, bp_mid, bp_to_high, P1_choice)

        # Otherwise, the array of mutations generated before is used as is
        else:
            mutations_nocost = mutations

        # The ancestral rates are copied ahead of the calculation of new mutations
        new_nocost = rates_nocost.copy()

        # Applying mutational effects depending on the epistasy that has been assumed, as before
        if epistasy == 'mult':
            new_nocost[:, 4:8] = new_nocost[:, 4:8] + (new_nocost[:, 4:8] * mutations_nocost)

        elif epistasy == 'add':
            new_nocost[:, 4:8] = new_nocost[:, 4:8] + (rates_init[:, 4:8] * mutations_nocost)

        # Computing mutant fitness
        fit_new_nocost = evol_funct.fit_noise_dupli(new_nocost[:, 4], new_nocost[:, 6], new_nocost[:, 5],
                                                    new_nocost[:, 7], new_nocost[:, 2], new_nocost[:, 1],
                                                    alpha_m, alpha_p, cv_0)

        # Ancestral and mutant fitness are normalized between 0 and 1
        fit_nocost = fit_nocost / 0.42
        fit_new_nocost = fit_new_nocost / 0.42

        # Calculating fixation probability from the fitness difference (ancestral - mutant)
        prob_nocost = fix_prob(fit_nocost, fit_new_nocost, pop_size)

        # Decisions on mutations
        nocost_decision = np.where(fix_decision < prob_nocost, 1, 0)

        # If the mutation takes Bm or Bp above or below the boundaries, it is canceled
        low_bm1 = np.where(new_nocost[:, 4] < bm_min, 0, 1)
        low_bm2 = np.where(new_nocost[:, 6] < bm_min, 0, 1)
        low_bp1 = np.where(new_nocost[:, 5] < bp_min, 0, 1)
        low_bp2 = np.where(new_nocost[:, 7] < bp_min, 0, 1)
        nocost_decision = nocost_decision * (low_bm1 * low_bm2 * low_bp1 * low_bp2)

        high_bm1 = np.where(new_nocost[:, 4] > bm_max, 0, 1)
        high_bm2 = np.where(new_nocost[:, 6] > bm_max, 0, 1)
        high_bp1 = np.where(new_nocost[:, 5] > bp_max, 0, 1)
        high_bp2 = np.where(new_nocost[:, 7] > bp_max, 0, 1)
        nocost_decision = nocost_decision * (high_bm1 * high_bm2 * high_bp1 * high_bp2)

        # Identifying all couples which are not mutated in the current round
        kept_nocost = np.abs(nocost_decision - 1)

        # Counting new mutations in copies P1 and P2, and adding them to the previous totals
        mut_P1_nocost = P1_choice * nocost_decision
        mut_P2_nocost = P2_choice * nocost_decision

        muta_nocost_current['P1 Mutations'] += mut_P1_nocost
        muta_nocost_current['P2 Mutations'] += mut_P2_nocost

        # Creating the array where only the mutants are kept
        mutants_nocost = new_nocost.copy()
        mutants_nocost[:, 4:8] = mutants_nocost[:, 4:8] * nocost_decision[:, np.newaxis]

        # Creating the array where only the non-mutated couples are kept
        nomut_nocost = rates_nocost.copy()
        nomut_nocost[:, 4:8] = nomut_nocost[:, 4:8] * kept_nocost[:, np.newaxis]

        # Combining them in a final array
        final_nocost = mutants_nocost.copy()
        final_nocost[:, 4:8] = final_nocost[:, 4:8] + nomut_nocost[:, 4:8]

    else:
        final_nocost = rates_nocost.copy()

        if dupli_loss:
            loss_data.at[0, 'No Cost'] = n_nocost

    # Third, for the ADS-only model
    if not done_ADS:

        if dupli_loss:
            # Gene loss is performed when neutral
            loss_ADS = evol_funct.gene_loss(rates_ADS, fit_ADS, evol_funct.fit_parabola, pop_size,
                                            (alpha_m, alpha_p))

            # The rates array is updated to take into account the newly performed gene loss
            rates_ADS = loss_ADS[0].copy()

            # The remaining number of intact couples is saved
            n_ADS = loss_ADS[1]
            loss_data.at[0, 'ADS'] = n_ADS

        # If bidirectional asymmetry of mutational effects has been used, the final array of mutations for this
        # simulation is generated
        if bidirectional:
            if bivariate:
                # Because this has not been implemented for a bivariate distribution of mutational effects
                raise Exception

            mutations_ADS = evol_funct.muta_bidirectional(rates_ADS, n_couples, (low_lim_bm, high_lim_bm),
                                                          (low_lim_bp, high_lim_bp), bm_to_low, bm_mid, bm_to_high,
                                                          bp_to_low, bp_mid, bp_to_high, P1_choice)

        # Otherwise, the array of mutations generated before is used as is
        else:
            mutations_ADS = mutations

        # The ancestral rates are copied ahead of the calculation of new mutations
        new_ADS = rates_ADS.copy()

        # Applying mutational effects
        if epistasy == 'mult':
            new_ADS[:, 4:8] = new_ADS[:, 4:8] + (new_ADS[:, 4:8] * mutations_ADS)

        elif epistasy == 'add':
            new_ADS[:, 4:8] = new_ADS[:, 4:8] + (start_ADS[:, 4:8] * mutations_ADS)

        # Computing mutant fitness
        fit_new_ADS = evol_funct.fit_parabola(new_ADS[:, 4], new_ADS[:, 6], new_ADS[:, 5], new_ADS[:, 7],
                                              new_ADS[:, 2], new_ADS[:, 1], alpha_m, alpha_p)

        # Ancestral and mutant fitness are normalized between 0 and 1
        fit_ADS = fit_ADS / 0.42
        fit_new_ADS = fit_new_ADS / 0.42

        # Calculating fixation probability from the fitness difference (ancestral - mutant)
        prob_ADS = fix_prob(fit_ADS, fit_new_ADS, pop_size)

        # Decisions on mutations
        ADS_decision = np.where(fix_decision < prob_ADS, 1, 0)

        # If the mutation takes Bm or Bp above or below the boundaries, it is canceled
        low_bm1 = np.where(new_ADS[:, 4] < bm_min, 0, 1)
        low_bm2 = np.where(new_ADS[:, 6] < bm_min, 0, 1)
        low_bp1 = np.where(new_ADS[:, 5] < bp_min, 0, 1)
        low_bp2 = np.where(new_ADS[:, 7] < bp_min, 0, 1)
        ADS_decision = ADS_decision * (low_bm1 * low_bm2 * low_bp1 * low_bp2)

        high_bm1 = np.where(new_ADS[:, 4] > bm_max, 0, 1)
        high_bm2 = np.where(new_ADS[:, 6] > bm_max, 0, 1)
        high_bp1 = np.where(new_ADS[:, 5] > bp_max, 0, 1)
        high_bp2 = np.where(new_ADS[:, 7] > bp_max, 0, 1)
        ADS_decision = ADS_decision * (high_bm1 * high_bm2 * high_bp1 * high_bp2)

        # Identifying all couples which are not mutated in the current round
        kept_ADS = np.abs(ADS_decision - 1)

        # Counting new mutations in copies P1 and P2, and adding them to the previous totals
        mut_P1_ADS = P1_choice * ADS_decision
        mut_P2_ADS = P2_choice * ADS_decision

        muta_ADS_current['P1 Mutations'] += mut_P1_ADS
        muta_ADS_current['P2 Mutations'] += mut_P2_ADS

        # Creating the array where only the mutants are kept
        mutants_ADS = new_ADS.copy()
        mutants_ADS[:, 4:8] = mutants_ADS[:, 4:8] * ADS_decision[:, np.newaxis]

        # Creating the array where only the non-mutated couples are kept
        nomut_ADS = rates_ADS.copy()
        nomut_ADS[:, 4:8] = nomut_ADS[:, 4:8] * kept_ADS[:, np.newaxis]

        # Combining them in a final array
        final_ADS = mutants_ADS.copy()
        final_ADS[:, 4:8] = final_ADS[:, 4:8] + nomut_ADS[:, 4:8]

        # Preparing a df to save
        to_save_ADS = final_ADS.copy()

    else:
        final_ADS = rates_ADS.copy()

        if dupli_loss:
            loss_data.at[0, 'ADS'] = n_ADS

    # Finally, for the "minimal" simulation
    if not done_min:

        if dupli_loss:
            # Gene loss is performed when neutral
            loss_min = evol_funct.gene_loss(rates_min, fit_min, evol_funct.fit_parabola, pop_size,
                                            (alpha_m, alpha_p))

            # The rates array is updated to take into account the newly performed gene loss
            rates_min = loss_min[0].copy()

            # The remaining number of intact couples is saved
            n_min = loss_min[1]
            loss_data.at[0, 'Minimal'] = n_min

        # If bidirectional asymmetry of mutational effects has been used, the final array of mutations for this
        # simulation is generated
        if bidirectional:
            if bivariate:
                # Because this has not been implemented for a bivariate distribution of mutational effects
                raise Exception

            mutations_min = evol_funct.muta_bidirectional(rates_min, n_couples, (low_lim_bm, high_lim_bm),
                                                          (low_lim_bp, high_lim_bp), bm_to_low, bm_mid, bm_to_high,
                                                          bp_to_low, bp_mid, bp_to_high, P1_choice)

        # Otherwise, the array of mutations generated before is used as is
        else:
            mutations_min = mutations

        # The ancestral rates are copied ahead of the calculation of new mutations
        new_min = rates_min.copy()

        # Applying mutational effects
        if epistasy == 'mult':
            new_min[:, 4:8] = new_min[:, 4:8] + (new_min[:, 4:8] * mutations_min)

        elif epistasy == 'add':
            new_min[:, 4:8] = new_min[:, 4:8] + (start_min[:, 4:8] * mutations_min)

        # Computing mutant fitness
        fit_new_min = evol_funct.fit_parabola(new_min[:, 4], new_min[:, 6], new_min[:, 5], new_min[:, 7],
                                              new_min[:, 2], new_min[:, 1], alpha_m, alpha_p)

        # Ancestral and mutant fitness are normalized between 0 and 1
        fit_min = fit_min / 0.42
        fit_new_min = fit_new_min / 0.42

        # Calculating fixation probability from the fitness difference (ancestral - mutant)
        prob_min = fix_prob(fit_min, fit_new_min, pop_size)

        # Decisions on mutations
        min_decision = np.where(fix_decision < prob_min, 1, 0)

        # If the mutation takes Bm or Bp above or below the boundaries, it is canceled
        low_bm1 = np.where(new_min[:, 4] < bm_min, 0, 1)
        low_bm2 = np.where(new_min[:, 6] < bm_min, 0, 1)
        low_bp1 = np.where(new_min[:, 5] < bp_min, 0, 1)
        low_bp2 = np.where(new_min[:, 7] < bp_min, 0, 1)
        min_decision = min_decision * (low_bm1 * low_bm2 * low_bp1 * low_bp2)

        high_bm1 = np.where(new_min[:, 4] > bm_max, 0, 1)
        high_bm2 = np.where(new_min[:, 6] > bm_max, 0, 1)
        high_bp1 = np.where(new_min[:, 5] > bp_max, 0, 1)
        high_bp2 = np.where(new_min[:, 7] > bp_max, 0, 1)
        min_decision = min_decision * (high_bm1 * high_bm2 * high_bp1 * high_bp2)

        # Identifying all couples which are not mutated in the current round
        kept_min = np.abs(min_decision - 1)

        # Counting new mutations in copies P1 and P2, and adding them to the previous totals
        mut_P1_min = P1_choice * min_decision
        mut_P2_min = P2_choice * min_decision

        muta_min_current['P1 Mutations'] += mut_P1_min
        muta_min_current['P2 Mutations'] += mut_P2_min

        # Creating the array where only the mutants are kept
        mutants_min = new_min.copy()
        mutants_min[:, 4:8] = mutants_min[:, 4:8] * min_decision[:, np.newaxis]

        # Creating the array where only the non-mutated couples are kept
        nomut_min = rates_min.copy()
        nomut_min[:, 4:8] = nomut_min[:, 4:8] * kept_min[:, np.newaxis]

        # Combining them in a final array
        final_min = mutants_min.copy()
        final_min[:, 4:8] = final_min[:, 4:8] + nomut_min[:, 4:8]

        # Preparing a df to save
        to_save_min = final_min.copy()

    else:
        final_min = rates_min.copy()

        if dupli_loss:
            loss_data.at[0, 'Minimal'] = n_min

    # The "final" couples are set as the new ancestral state, in prevision of the next mutation round
    rates_mixed = final_mixed.copy()
    rates_nocost = final_nocost.copy()
    rates_ADS = final_ADS.copy()
    rates_min = final_min.copy()

    # The data on the number of intact couples is saved (if appropriate)
    if dupli_loss:
        loss_data.to_csv(os.path.join(path_folder, f'Gene_loss_{run_name}.csv'), index=False, mode='a', header=False)

    # Dataframe of expression rates and protein abundance are made for the four selection regimes
    data_model['Round'] = step

    if not done_mixed:
        data_mixed = data_model.copy()
        data_mixed.iloc[:, 2:9] = rates_mixed[:, 1:8]
        data_mixed['Prot1'] = (data_mixed['Bm1'] * data_mixed['Bp1']) / (alpha_m * alpha_p)
        data_mixed['Prot2'] = (data_mixed['Bm2'] * data_mixed['Bp2']) / (alpha_m * alpha_p)
        data_mixed['cv1'] = np.sqrt((1/data_mixed['Prot1']) + (alpha_p/data_mixed['Bm1']) + cv_0**2)
        data_mixed['cv2'] = np.sqrt((1/data_mixed['Prot2']) + (alpha_p/data_mixed['Bm2']) + cv_0**2)
        data_mixed['Exp_cost'] = lm * c_m * (data_mixed['Bm1'] + data_mixed['Bm2'])
        data_mixed['Model'] = 'Mixed'

    if not done_nocost:
        data_nocost = data_model.copy()
        data_nocost.iloc[:, 2:9] = rates_nocost[:, 1:8]
        data_nocost['Prot1'] = (data_nocost['Bm1'] * data_nocost['Bp1']) / (alpha_m * alpha_p)
        data_nocost['Prot2'] = (data_nocost['Bm2'] * data_nocost['Bp2']) / (alpha_m * alpha_p)
        data_nocost['cv1'] = np.sqrt((1/data_nocost['Prot1']) + (alpha_p/data_nocost['Bm1']) + cv_0**2)
        data_nocost['cv2'] = np.sqrt((1/data_nocost['Prot2']) + (alpha_p/data_nocost['Bm2']) + cv_0**2)
        data_nocost['Exp_cost'] = lm * c_m * (data_nocost['Bm1'] + data_nocost['Bm2'])
        data_nocost['Model'] = 'No Cost'

    if not done_ADS:
        data_ADS = data_model.copy()
        data_ADS.iloc[:, 2:9] = rates_ADS[:, 1:8]
        data_ADS['Prot1'] = (data_ADS['Bm1'] * data_ADS['Bp1']) / (alpha_m * alpha_p)
        data_ADS['Prot2'] = (data_ADS['Bm2'] * data_ADS['Bp2']) / (alpha_m * alpha_p)
        data_ADS['cv1'] = np.sqrt((1/data_ADS['Prot1']) + (alpha_p/data_ADS['Bm1']) + cv_0**2)
        data_ADS['cv2'] = np.sqrt((1/data_ADS['Prot2']) + (alpha_p/data_ADS['Bm2']) + cv_0**2)
        data_ADS['Exp_cost'] = lm * c_m * (data_ADS['Bm1'] + data_ADS['Bm2'])
        data_ADS['Model'] = 'ADS'

    if not done_min:
        data_min = data_model.copy()
        data_min.iloc[:, 2:9] = rates_min[:, 1:8]
        data_min['Prot1'] = (data_min['Bm1'] * data_min['Bp1']) / (alpha_m * alpha_p)
        data_min['Prot2'] = (data_min['Bm2'] * data_min['Bp2']) / (alpha_m * alpha_p)
        data_min['cv1'] = np.sqrt((1/data_min['Prot1']) + (alpha_p/data_min['Bm1']) + cv_0**2)
        data_min['cv2'] = np.sqrt((1/data_min['Prot2']) + (alpha_p/data_min['Bm2']) + cv_0**2)
        data_min['Exp_cost'] = lm * c_m * (data_min['Bm1'] + data_min['Bm2'])
        data_min['Model'] = 'Minimal'

    # The level of divergence between duplicates is assessed.
    # The equality of the medians of logfold divergence in protein abundance is tested.
    p_values = []
    for regime in [data_mixed, data_nocost, data_ADS, data_min]:
        fold_log = evol_funct.fold_change('Prot1', 'Prot2', regime)
        value = stats.median_test(fold_log, fold_med[pEst_fold_param], nan_policy='omit')[1]

        p_values.append(value)

    pval_current = Mood_pvals_model.copy()
    pval_current['Round'] = step
    pval_current.iloc[0, 1:5] = p_values

    Mood_pvals = pd.concat([Mood_pvals, pval_current]).reset_index(drop=True)

    #pval_current.to_csv(os.path.join(path_folder, f'Mood_pvalues_{run_name}.csv'), mode='a', header=False, index=False)

    # The number of fixed mutations per paralog copy is saved for the four models for some pre-specified rounds.
    # This is only done while the end condition is not true
    if step in [10, 50, 100] or step % 500 == 0:
        muta_mixed_current['Round'] = step
        muta_nocost_current['Round'] = step
        muta_ADS_current['Round'] = step
        muta_min_current['Round'] = step

        if not done_mixed:
            muta_mixed = pd.concat([muta_mixed, muta_mixed_current]).reset_index(drop=True)

        if not done_nocost:
            muta_nocost = pd.concat([muta_nocost, muta_nocost_current]).reset_index(drop=True)

        if not done_ADS:
            muta_ADS = pd.concat([muta_ADS, muta_ADS_current]).reset_index(drop=True)

        if not done_min:
            muta_min = pd.concat([muta_min, muta_min_current]).reset_index(drop=True)

    # The end conditions are evaluated. If it is the first time that they are true, the data for the current round
    # (dataframes made previously) is saved. The number of fixed mutations is also saved again if the end condition is
    # verified for the first time.
    if end_cond == 'pEst':

        # For the mixed model
        if p_values[0] >= 0.1:
            last_mixed = True

        # For the no-cost model
        if p_values[1] >= 0.1:
            last_nocost = True

        # For the ADS-only model
        if p_values[2] >= 0.1:
            last_ADS = True

        # For the "minimal" model
        if p_values[3] >= 0.1:
            last_min = True

    elif end_cond == 'loss':

        # For the mixed model
        if n_mixed <= 0.15 * n_couples:
            last_mixed = True

        # For the no-cost model
        if n_nocost <= 0.15 * n_couples:
            last_nocost = True

        # For the ADS-only model
        if n_ADS <= 0.15 * n_couples:
            last_ADS = True

        # For the ADS-only model
        if n_min <= 0.15 * n_couples:
            last_min = True

    # If full-data saving has been enabled, transcription and translation rates for the current round are saved
    if full_data:
        if not last_mixed and not done_mixed:
            data_mixed['Fitness'] = evol_funct.fit_global_dupli(data_mixed['Bm1'], data_mixed['Bm2'], data_mixed['Bp1'],
                                                                data_mixed['Bp2'], data_mixed['pOpt'], data_mixed['Q'],
                                                                alpha_m, alpha_p, cv_0, data_mixed['lm'], c_m)

            data_all = pd.concat([data_all, data_mixed]).reset_index(drop=True)
            #data_mixed.to_csv(os.path.join(path_folder, f'data_Mixed_{run_name}.csv'), mode='a', header=False,
                              #index=False)

        if not last_nocost and not done_nocost:
            data_nocost['Fitness'] = evol_funct.fit_noise_dupli(data_nocost['Bm1'], data_nocost['Bm2'],
                                                                data_nocost['Bp1'],
                                                                data_nocost['Bp2'], data_nocost['pOpt'],
                                                                data_nocost['Q'], alpha_m, alpha_p, cv_0)
            data_all = pd.concat([data_all, data_nocost]).reset_index(drop=True)
            #data_nocost.to_csv(os.path.join(path_folder, f'data_NoCost_{run_name}.csv'), mode='a', header=False,
                               #index=False)

        if not last_ADS and not done_ADS:
            data_ADS['Fitness'] = evol_funct.fit_parabola(data_ADS['Bm1'], data_ADS['Bm2'], data_ADS['Bp1'],
                                                          data_ADS['Bp2'], data_ADS['pOpt'], data_ADS['Q'],
                                                          alpha_m, alpha_p)
            data_all = pd.concat([data_all, data_ADS]).reset_index(drop=True)
            #data_ADS.to_csv(os.path.join(path_folder, f'data_ADS_{run_name}.csv'), mode='a', header=False, index=False)

        if not last_min and not done_min:
            data_min['Fitness'] = evol_funct.fit_parabola(data_min['Bm1'], data_min['Bm2'], data_min['Bp1'],
                                                          data_min['Bp2'], data_min['pOpt'], data_min['Q'],
                                                          alpha_m, alpha_p)
            data_all = pd.concat([data_all, data_min]).reset_index(drop=True)
            #data_min.to_csv(os.path.join(path_folder, f'data_minimal_{run_name}.csv'), mode='a', header=False,
                            #index=False)

    # Data saving if the end condition has newly been verified
    if last_mixed and not done_mixed:
        data_mixed['Fitness'] = evol_funct.fit_global_dupli(data_mixed['Bm1'], data_mixed['Bm2'], data_mixed['Bp1'],
                                                            data_mixed['Bp2'], data_mixed['pOpt'], data_mixed['Q'],
                                                            alpha_m, alpha_p, cv_0, data_mixed['lm'], c_m)
        data_all = pd.concat([data_all, data_mixed]).reset_index(drop=True)
        #data_mixed.to_csv(os.path.join(path_folder, f'data_Mixed_{run_name}.csv'), mode='a', header=False, index=False)

        muta_mixed_current['Round'] = step
        muta_mixed = pd.concat([muta_mixed, muta_mixed_current]).reset_index(drop=True)

        #muta_mixed.to_csv(os.path.join(path_folder, f'muta_mixed_{run_name}.csv'), index=False)

    if last_nocost and not done_nocost:
        data_nocost['Fitness'] = evol_funct.fit_noise_dupli(data_nocost['Bm1'], data_nocost['Bm2'],
                                                            data_nocost['Bp1'], data_nocost['Bp2'],
                                                            data_nocost['pOpt'], data_nocost['Q'], alpha_m,
                                                            alpha_p, cv_0)
        data_all = pd.concat([data_all, data_nocost]).reset_index(drop=True)
        #data_nocost.to_csv(os.path.join(path_folder, f'data_NoCost_{run_name}.csv'), header=False, index=False)

        muta_nocost_current['Round'] = step
        muta_nocost = pd.concat([muta_nocost, muta_nocost_current]).reset_index(drop=True)

        #muta_nocost.to_csv(os.path.join(path_folder, f'muta_nocost_{run_name}.csv'), index=False)

    if last_ADS and not done_ADS:
        data_ADS['Fitness'] = evol_funct.fit_parabola(data_ADS['Bm1'], data_ADS['Bm2'], data_ADS['Bp1'],
                                                      data_ADS['Bp2'], data_ADS['pOpt'], data_ADS['Q'],
                                                      alpha_m, alpha_p)
        data_all = pd.concat([data_all, data_ADS]).reset_index(drop=True)
        #data_ADS.to_csv(os.path.join(path_folder, f'data_ADS_{run_name}.csv'), header=False, index=False)

        muta_ADS_current['Round'] = step
        muta_ADS = pd.concat([muta_ADS, muta_ADS_current]).reset_index(drop=True)

        #muta_ADS.to_csv(os.path.join(path_folder, f'muta_ADS_{run_name}.csv'), index=False)

    if last_min and not done_min:
        data_min['Fitness'] = evol_funct.fit_parabola(data_min['Bm1'], data_min['Bm2'], data_min['Bp1'],
                                                      data_min['Bp2'], data_min['pOpt'], data_min['Q'],
                                                      alpha_m, alpha_p)
        data_all = pd.concat([data_all, data_min]).reset_index(drop=True)
        #data_min.to_csv(os.path.join(path_folder, f'data_minimal_{run_name}.csv'), mode='a', header=False, index=False)

        muta_min_current['Round'] = step
        muta_min = pd.concat([muta_min, muta_min_current]).reset_index(drop=True)

        #muta_min.to_csv(os.path.join(path_folder, f'muta_min_{run_name}.csv'), index=False)

    # The end conditions are officially evaluated
    if last_mixed:
        done_mixed = True

    if last_nocost:
        done_nocost = True

    if last_ADS:
        done_ADS = True

    if last_min:
        done_min = True

    # Global end condition
    if (done_mixed and done_nocost and done_ADS and done_min):
        stop_sim = True

    if step % 100 == 0:
        print(f'Done with step {step}')

    # Saving the dataframes once the simulation is finished
    if stop_sim:
        data_all.to_csv(os.path.join(path_folder, f'data_all_{run_name}.csv'), index=False)

        muta_mixed['Model'] = 'Mixed'
        muta_nocost['Model'] = 'No Cost'
        muta_ADS['Model'] = 'ADS'
        muta_min['Model'] = 'Minimal'

        muta_all = pd.concat([muta_mixed, muta_nocost, muta_ADS, muta_min]).reset_index(drop=True)
        muta_all.to_csv(os.path.join(path_folder, f'muta_all_{run_name}.csv'), index=False)

        Mood_pvals.to_csv(os.path.join(path_folder, f'Mood_pvalues_{run_name}.csv'), index=False)

    step += 1

# End time for the mutation-selection rounds only
t1_evol = time.process_time()

# Once the simulation is completed, new datasets are first generated for each model (the two main models and the two
# internal controls). They contain only the "true duplicate" pairs, that is the subset of simulated pairs for which
# the loss of a gene copy would still not be tolerated by selection at the end of the simulation
dict_mixed = {'alpha_m': alpha_m, 'alpha_p': alpha_p, 'cv_0': cv_0, 'lm': rates_mixed[:, 3], 'c_m': c_m}
keys_mixed = ['alpha_m', 'alpha_p', 'cv_0', 'lm', 'c_m']
data_true_mixed = evol_funct.extract_true_dupli(rates_mixed, data_mixed, evol_funct.fit_global_dupli, dict_mixed,
                                                keys_mixed, pop_size, data_model, 'Mixed', f'data_all_{run_name}.csv',
                                                path_folder)

dict_nocost = {'alpha_m': alpha_m, 'alpha_p': alpha_p, 'cv_0': cv_0, 'lm': rates_nocost[:, 3], 'c_m': c_m}
keys_nocost = ['alpha_m', 'alpha_p', 'cv_0']
data_true_nocost = evol_funct.extract_true_dupli(rates_nocost, data_nocost, evol_funct.fit_noise_dupli, dict_nocost,
                                                 keys_nocost, pop_size, data_model, 'No Cost',
                                                 f'data_all_{run_name}.csv', path_folder)

dict_ADS = {'alpha_m': alpha_m, 'alpha_p': alpha_p, 'cv_0': cv_0, 'lm': rates_ADS[:, 3], 'c_m': c_m}
keys_ADS = ['alpha_m', 'alpha_p']
data_true_ADS = evol_funct.extract_true_dupli(rates_ADS, data_ADS, evol_funct.fit_parabola, dict_ADS,
                                              keys_ADS, pop_size, data_model, 'ADS', f'data_all_{run_name}.csv',
                                              path_folder)

dict_min = {'alpha_m': alpha_m, 'alpha_p': alpha_p, 'cv_0': cv_0, 'lm': rates_min[:, 3], 'c_m': c_m}
keys_min = ['alpha_m', 'alpha_p']
data_true_min = evol_funct.extract_true_dupli(rates_min, data_min, evol_funct.fit_parabola, dict_min,
                                              keys_min, pop_size, data_model, 'Minimal', f'data_all_{run_name}.csv',
                                              path_folder)

data_true_all = pd.concat([data_true_mixed, data_true_nocost, data_true_ADS, data_true_min]).reset_index(drop=True)
data_true_all.to_csv(os.path.join(path_folder, f'data_true_all_{run_name}.csv'), index=False)

# Log2 fold-changes are then computed for the two datasets (all pairs and true duplicates) corresponding to each of the
# four models.

# Final data is selected
end_mixed = data_mixed.copy()
end_nocost = data_nocost.copy()
end_ADS = data_ADS.copy()
end_min = data_min.copy()

end_true_mixed = evol_funct.assign_paralogs(data_true_mixed)[2]
end_true_nocost = evol_funct.assign_paralogs(data_true_nocost)[2]
end_true_ADS = evol_funct.assign_paralogs(data_true_ADS)[2]
end_true_min = evol_funct.assign_paralogs(data_true_min)[2]

end_true_mixed = end_true_mixed.reset_index(drop=True)
end_true_nocost = end_true_nocost.reset_index(drop=True)
end_true_ADS = end_true_ADS.reset_index(drop=True)
end_true_min = end_true_min.reset_index(drop=True)

# Log2 fold changes between duplicates are calculated for transcription rate, translation rate and protein abundance;
# Prior to this, singletons are removed. For the "true duplicates" set, this has already been done, but zeros
# are dropped again to prevent any problems with the calculation of log2 fold-changes:
for end_df in [end_mixed, end_nocost, end_ADS, end_min,
               end_true_mixed, end_true_nocost, end_true_ADS, end_true_min]:

    end_df.iloc[:, 5:] = end_df.iloc[:, 5:].mask(end_df.iloc[:, 5:] == 0)  # If Bm1==0
    end_df.iloc[:, 7:] = end_df.iloc[:, 7:].mask(end_df.iloc[:, 7:] == 0)  # If Bm2==0
    end_df.iloc[:, 6:] = end_df.iloc[:, 6:].mask(end_df.iloc[:, 6:] == 0)  # If Bp1==0
    end_df.iloc[:, 8:] = end_df.iloc[:, 8:].mask(end_df.iloc[:, 8:] == 0)  # If Bp2==0

end_mixed = end_mixed.dropna(axis=0, subset=['Bm1', 'Bp1', 'Bm2', 'Bp2']).reset_index(drop=True)
end_true_mixed = end_true_mixed.dropna(axis=0, subset=['Bm1', 'Bp1', 'Bm2', 'Bp2']).reset_index(drop=True)

end_nocost = end_nocost.dropna(axis=0, subset=['Bm1', 'Bp1', 'Bm2', 'Bp2']).reset_index(drop=True)
end_true_nocost = end_true_nocost.dropna(axis=0, subset=['Bm1', 'Bp1', 'Bm2', 'Bp2']).reset_index(drop=True)

end_ADS = end_ADS.dropna(axis=0, subset=['Bm1', 'Bp1', 'Bm2', 'Bp2']).reset_index(drop=True)
end_true_ADS = end_true_ADS.dropna(axis=0, subset=['Bm1', 'Bp1', 'Bm2', 'Bp2']).reset_index(drop=True)

end_min = end_min.dropna(axis=0, subset=['Bm1', 'Bp1', 'Bm2', 'Bp2']).reset_index(drop=True)
end_true_min = end_true_min.dropna(axis=0, subset=['Bm1', 'Bp1', 'Bm2', 'Bp2']).reset_index(drop=True)

# Then, the log2 fold-changes are computed, as well as the divergence ratios
for end_df in [end_mixed, end_nocost, end_ADS, end_min, end_true_mixed, end_true_nocost, end_true_ADS, end_true_min]:

    logfold_bm = evol_funct.fold_change('Bm1', 'Bm2', end_df)
    logfold_bp = evol_funct.fold_change('Bp1', 'Bp2', end_df)
    logfold_prot = evol_funct.fold_change('Prot1', 'Prot2', end_df)
    div_ratio = np.log2(evol_funct.div_ratio('Bm1', 'Bm2', 'Bp1', 'Bp2', end_df)['Divergence_ratio'])

    end_df.insert(16, 'Transcription rate', logfold_bm)
    end_df.insert(17, 'Translation rate', logfold_bp)
    end_df.insert(18, 'Protein abundance', logfold_prot)
    end_df.insert(19, 'Divergence ratio', div_ratio)

end_mixed_df = end_mixed.iloc[:, 15:20].copy()
end_nocost_df = end_nocost.iloc[:, 15:20].copy()
end_ADS_df = end_ADS.iloc[:, 15:20].copy()
end_min_df = end_min.iloc[:, 15:20].copy()

end_true_mixed_df = end_true_mixed.iloc[:, 15:20].copy()
end_true_nocost_df = end_true_nocost.iloc[:, 15:20].copy()
end_true_ADS_df = end_true_ADS.iloc[:, 15:20].copy()
end_true_min_df = end_true_min.iloc[:, 15:20].copy()

fold_yeast = folds_df[['bm_fold_original', bp_fold_param, pEst_fold_param]].copy()
fold_yeast.insert(0, 'Model', 'Yeast Data')
fold_yeast.columns = ['Model', 'Transcription rate', 'Translation rate', 'Protein abundance']

div_yeast = np.log2(evol_funct.div_ratio_log('bm_P1', 'bm_P2', f'{bp_param}_P1', f'{bp_param}_P2',
                                             folds_df)['Divergence_ratio'])
fold_yeast.insert(4, 'Divergence ratio', div_yeast)

# The data for WGD or SSD-derived yeast couples is handled similarly as the full yeast dataset
couples_WGD = folds_df[folds_df['Duplication'] == 'WGD'][['bm_fold_original', bp_fold_param, pEst_fold_param]].copy()
couples_WGD.insert(0, 'Model', 'WGD couples')
couples_WGD.columns = ['Model', 'Transcription rate', 'Translation rate', 'Protein abundance']
couples_WGD = couples_WGD.reset_index(drop=True)

folds_WGD = folds_df[folds_df['Duplication'] == 'WGD'].reset_index(drop=True)
div_WGD = np.log2(evol_funct.div_ratio_log('bm_P1', 'bm_P2', f'{bp_param}_P1', f'{bp_param}_P2',
                                           folds_WGD)['Divergence_ratio'])
couples_WGD.insert(4, 'Divergence ratio', div_WGD)

couples_SSD = folds_df[folds_df['Duplication'] == 'SSD'][['bm_fold_original', bp_fold_param, pEst_fold_param]].copy()
couples_SSD.insert(0, 'Model', 'SSD couples')
couples_SSD.columns = ['Model', 'Transcription rate', 'Translation rate', 'Protein abundance']
couples_SSD = couples_SSD.reset_index(drop=True)

folds_SSD = folds_df[folds_df['Duplication'] == 'SSD'].reset_index(drop=True)
div_SSD = np.log2(evol_funct.div_ratio_log('bm_P1', 'bm_P2', f'{bp_param}_P1', f'{bp_param}_P2',
                                           folds_SSD)['Divergence_ratio'])
couples_SSD.insert(4, 'Divergence ratio', div_SSD)

# Two dataframes (one for all simulated pairs and the other only for true duplicates) are made combining this data.
# These will mostly be useful for the later construction of figures
folds_all = pd.concat([end_mixed_df, end_nocost_df, end_ADS_df, end_min_df, fold_yeast, couples_WGD, couples_SSD]).reset_index(drop=True)
folds_melted = folds_all.melt(id_vars=['Model'], value_vars=['Transcription rate', 'Translation rate',
                                                             'Protein abundance', 'Divergence ratio'],
                              var_name='Parameter', value_name='Log2 fold-change')
folds_melted.to_csv(os.path.join(path_folder, f'Fold_changes_all_{run_name}.csv'), index=False)

folds_true = pd.concat([end_true_mixed_df, end_true_nocost_df, end_true_ADS_df, end_true_min_df,
                        fold_yeast, couples_WGD, couples_SSD])
folds_melt_true = folds_true.melt(id_vars=['Model'], value_vars=['Transcription rate', 'Translation rate',
                                                                 'Protein abundance', 'Divergence ratio'],
                                  var_name='Parameter', value_name='Log2 fold-change')
folds_melt_true.to_csv(os.path.join(path_folder, f'Fold_changes_true_{run_name}.csv'), index=False)

# KS tests as well as Mood's median tests are performed to compare all divergence distributions. The comparisons are
# made between each simulated distribution and the yeast data. The results from these calculations are combined into
# a single dataframe, which is then exported as a csv
dist_model = pd.DataFrame(columns=['Model', 'Mut_sigma', 'Mut_ratio', 'Mut_alpha', 'Mut_corr', 'Iter',
                                  'Folder', 'Set', 'Comparison', 'Property', 'Type', 'Moods_stat', 'Moods_p-val',
                                  'KS_stat', 'KS_p-val'])

dist_comp = dist_model.copy()

row_num = 0
for model in ['Mixed', 'No Cost', 'ADS', 'Minimal']:
    for dupli_prop in ['Transcription rate', 'Translation rate', 'Protein abundance', 'Divergence ratio']:

        # 1) Comparisons with all yeast duplicates
        # First for all pairs
        folds_subset = folds_all[folds_all['Model'] == model]
        Moods_test = stats.median_test(folds_subset[dupli_prop], fold_yeast[dupli_prop])
        KS_test = stats.ks_2samp(folds_subset[dupli_prop], fold_yeast[dupli_prop])

        all_all = [model, float(bm_params[1]), ratio, skew_param, correlation, run_name, sim_folder, sim_set,
                   'All duplicates', dupli_prop, 'All', Moods_test[0], Moods_test[1], KS_test[0], KS_test[1]]
        dist_current = dist_model.copy()
        dist_current.loc[0] = all_all
        dist_comp = pd.concat([dist_comp, dist_current]).reset_index(drop=True)

        # Then for true duplicates
        folds_subset = folds_true[folds_true['Model'] == model]
        Moods_test = stats.median_test(folds_subset[dupli_prop], fold_yeast[dupli_prop])
        KS_test = stats.ks_2samp(folds_subset[dupli_prop], fold_yeast[dupli_prop])

        all_true = [model, float(bm_params[1]), ratio, skew_param, correlation, run_name, sim_folder, sim_set,
                   'All duplicates', dupli_prop, 'True duplicates', Moods_test[0], Moods_test[1], KS_test[0], KS_test[1]]
        dist_current = dist_model.copy()
        dist_current.loc[0] = all_true
        dist_comp = pd.concat([dist_comp, dist_current]).reset_index(drop=True)

        # 2) Comparisons with WGD-derived paralogs
        # First for all pairs
        folds_subset = folds_all[folds_all['Model'] == model]
        Moods_test = stats.median_test(folds_subset[dupli_prop], couples_WGD[dupli_prop])
        KS_test = stats.ks_2samp(folds_subset[dupli_prop], couples_WGD[dupli_prop])

        all_all = [model, float(bm_params[1]), ratio, skew_param, correlation, run_name, sim_folder, sim_set,
                   'WGD', dupli_prop, 'All', Moods_test[0], Moods_test[1], KS_test[0], KS_test[1]]
        dist_current = dist_model.copy()
        dist_current.loc[0] = all_all
        dist_comp = pd.concat([dist_comp, dist_current]).reset_index(drop=True)

        # Then for true duplicates
        folds_subset = folds_true[folds_true['Model'] == model]
        Moods_test = stats.median_test(folds_subset[dupli_prop], couples_WGD[dupli_prop])
        KS_test = stats.ks_2samp(folds_subset[dupli_prop], couples_WGD[dupli_prop])

        all_true = [model, float(bm_params[1]), ratio, skew_param, correlation, run_name, sim_folder, sim_set,
                    'WGD', dupli_prop, 'True duplicates', Moods_test[0], Moods_test[1], KS_test[0], KS_test[1]]
        dist_current = dist_model.copy()
        dist_current.loc[0] = all_true
        dist_comp = pd.concat([dist_comp, dist_current]).reset_index(drop=True)

        # 3) Comparisons with SSD-derived paralogs
        # First for all simulated pairs
        folds_subset = folds_all[folds_all['Model'] == model]
        Moods_test = stats.median_test(folds_subset[dupli_prop], couples_SSD[dupli_prop])
        KS_test = stats.ks_2samp(folds_subset[dupli_prop], couples_SSD[dupli_prop])

        all_all = [model, float(bm_params[1]), ratio, skew_param, correlation, run_name, sim_folder, sim_set,
                   'SSD', dupli_prop, 'All', Moods_test[0], Moods_test[1], KS_test[0], KS_test[1]]
        dist_current = dist_model.copy()
        dist_current.loc[0] = all_all
        dist_comp = pd.concat([dist_comp, dist_current]).reset_index(drop=True)

        # Then for true duplicates
        folds_subset = folds_true[folds_true['Model'] == model]
        Moods_test = stats.median_test(folds_subset[dupli_prop], couples_SSD[dupli_prop])
        KS_test = stats.ks_2samp(folds_subset[dupli_prop], couples_SSD[dupli_prop])

        all_true = [model, float(bm_params[1]), ratio, skew_param, correlation, run_name, sim_folder, sim_set,
                    'SSD', dupli_prop, 'True duplicates', Moods_test[0], Moods_test[1], KS_test[0], KS_test[1]]
        dist_current = dist_model.copy()
        dist_current.loc[0] = all_true
        dist_comp = pd.concat([dist_comp, dist_current]).reset_index(drop=True)

# Exporting the resulting dataframe
dist_comp.to_csv(os.path.join(path_folder, f'distribution_comps_{run_name}.csv'), index=False)

# Finally, the three relevant within-pair divergence correlations (Bm and Bp log-fold changes, Signed Bm and Bp log-fold
# changes and divergence ratio with log-fold protein abundance change) are computed, both for all pairs as well as
# for the subset of true duplicates.

# This first requires computing the signed fold-changes. They are computed in the two possible orientations, meaning
# that the dataset is effectively duplicated
folds_signed_model = pd.DataFrame(columns=['Model', 'Bm1', 'Bm2', 'Bp1', 'Bp2', 'Signed_Bm', 'Signed_Bp'])

# This is first done for the complete set of simulated gene pairs
# For the mixed model
mixed_sign_P1 = folds_signed_model.copy()
mixed_sign_P1[['Bm1', 'Bm2', 'Bp1', 'Bp2']] = end_mixed[['Bm1', 'Bm2', 'Bp1', 'Bp2']].copy()
mixed_sign_P1['Model'] = 'Mixed'
mixed_sign_P2 = mixed_sign_P1.copy()
mixed_sign_P1['Signed_Bm'] = np.log2(mixed_sign_P1['Bm1'] / mixed_sign_P1['Bm2'])
mixed_sign_P1['Signed_Bp'] = np.log2(mixed_sign_P1['Bp1'] / mixed_sign_P1['Bp2'])
mixed_sign_P2['Signed_Bm'] = np.log2(mixed_sign_P2['Bm2'] / mixed_sign_P2['Bm1'])
mixed_sign_P2['Signed_Bp'] = np.log2(mixed_sign_P2['Bp1'] / mixed_sign_P2['Bp1'])
mixed_signed = pd.concat([mixed_sign_P1, mixed_sign_P2]).reset_index(drop=True)

# For the no-cost simulation
nocost_sign_P1 = folds_signed_model.copy()
nocost_sign_P1[['Bm1', 'Bm2', 'Bp1', 'Bp2']] = end_nocost[['Bm1', 'Bm2', 'Bp1', 'Bp2']].copy()
nocost_sign_P1['Model'] = 'No Cost'
nocost_sign_P2 = nocost_sign_P1.copy()
nocost_sign_P1['Signed_Bm'] = np.log2(nocost_sign_P1['Bm1'] / nocost_sign_P1['Bm2'])
nocost_sign_P1['Signed_Bp'] = np.log2(nocost_sign_P1['Bp1'] / nocost_sign_P1['Bp2'])
nocost_sign_P2['Signed_Bm'] = np.log2(nocost_sign_P2['Bm2'] / nocost_sign_P2['Bm1'])
nocost_sign_P2['Signed_Bp'] = np.log2(nocost_sign_P2['Bp1'] / nocost_sign_P2['Bp1'])
nocost_signed = pd.concat([nocost_sign_P1, nocost_sign_P2]).reset_index(drop=True)

# For the ADS simulation
ADS_sign_P1 = folds_signed_model.copy()
ADS_sign_P1[['Bm1', 'Bm2', 'Bp1', 'Bp2']] = end_ADS[['Bm1', 'Bm2', 'Bp1', 'Bp2']].copy()
ADS_sign_P1['Model'] = 'ADS'
ADS_sign_P2 = ADS_sign_P1.copy()
ADS_sign_P1['Signed_Bm'] = np.log2(ADS_sign_P1['Bm1'] / ADS_sign_P1['Bm2'])
ADS_sign_P1['Signed_Bp'] = np.log2(ADS_sign_P1['Bp1'] / ADS_sign_P1['Bp2'])
ADS_sign_P2['Signed_Bm'] = np.log2(ADS_sign_P2['Bm2'] / ADS_sign_P2['Bm1'])
ADS_sign_P2['Signed_Bp'] = np.log2(ADS_sign_P2['Bp1'] / ADS_sign_P2['Bp1'])
ADS_signed = pd.concat([ADS_sign_P1, ADS_sign_P2]).reset_index(drop=True)

# For the minimal model
min_sign_P1 = folds_signed_model.copy()
min_sign_P1[['Bm1', 'Bm2', 'Bp1', 'Bp2']] = end_min[['Bm1', 'Bm2', 'Bp1', 'Bp2']].copy()
min_sign_P1['Model'] = 'Minimal'
min_sign_P2 = min_sign_P1.copy()
min_sign_P1['Signed_Bm'] = np.log2(min_sign_P1['Bm1'] / min_sign_P1['Bm2'])
min_sign_P1['Signed_Bp'] = np.log2(min_sign_P1['Bp1'] / min_sign_P1['Bp2'])
min_sign_P2['Signed_Bm'] = np.log2(min_sign_P2['Bm2'] / min_sign_P2['Bm1'])
min_sign_P2['Signed_Bp'] = np.log2(min_sign_P2['Bp1'] / min_sign_P2['Bp1'])
min_signed = pd.concat([min_sign_P1, min_sign_P2]).reset_index(drop=True)

# Then, the same is done when considering only the set of true duplicates:
# For the mixed model
mixed_true_sign_P1 = folds_signed_model.copy()
mixed_true_sign_P1[['Bm1', 'Bm2', 'Bp1', 'Bp2']] = end_true_mixed[['Bm1', 'Bm2', 'Bp1', 'Bp2']].copy()
mixed_true_sign_P1['Model'] = 'Mixed'
mixed_true_sign_P2 = mixed_sign_P1.copy()
mixed_true_sign_P1['Signed_Bm'] = np.log2(mixed_true_sign_P1['Bm1'] / mixed_true_sign_P1['Bm2'])
mixed_true_sign_P1['Signed_Bp'] = np.log2(mixed_true_sign_P1['Bp1'] / mixed_true_sign_P1['Bp2'])
mixed_true_sign_P2['Signed_Bm'] = np.log2(mixed_true_sign_P2['Bm2'] / mixed_true_sign_P2['Bm1'])
mixed_true_sign_P2['Signed_Bp'] = np.log2(mixed_true_sign_P2['Bp1'] / mixed_true_sign_P2['Bp1'])
mixed_true_signed = pd.concat([mixed_true_sign_P1, mixed_true_sign_P2]).reset_index(drop=True)

# For the no-cost simulation
nocost_true_sign_P1 = folds_signed_model.copy()
nocost_true_sign_P1[['Bm1', 'Bm2', 'Bp1', 'Bp2']] = end_true_nocost[['Bm1', 'Bm2', 'Bp1', 'Bp2']].copy()
nocost_true_sign_P1['Model'] = 'No Cost'
nocost_true_sign_P2 = nocost_sign_P1.copy()
nocost_true_sign_P1['Signed_Bm'] = np.log2(nocost_true_sign_P1['Bm1'] / nocost_true_sign_P1['Bm2'])
nocost_true_sign_P1['Signed_Bp'] = np.log2(nocost_true_sign_P1['Bp1'] / nocost_true_sign_P1['Bp2'])
nocost_true_sign_P2['Signed_Bm'] = np.log2(nocost_true_sign_P2['Bm2'] / nocost_true_sign_P2['Bm1'])
nocost_true_sign_P2['Signed_Bp'] = np.log2(nocost_true_sign_P2['Bp1'] / nocost_true_sign_P2['Bp1'])
nocost_true_signed = pd.concat([nocost_true_sign_P1, nocost_true_sign_P2]).reset_index(drop=True)

# For the ADS simulation
ADS_true_sign_P1 = folds_signed_model.copy()
ADS_true_sign_P1[['Bm1', 'Bm2', 'Bp1', 'Bp2']] = end_true_ADS[['Bm1', 'Bm2', 'Bp1', 'Bp2']].copy()
ADS_true_sign_P1['Model'] = 'ADS'
ADS_true_sign_P2 = ADS_sign_P1.copy()
ADS_true_sign_P1['Signed_Bm'] = np.log2(ADS_true_sign_P1['Bm1'] / ADS_true_sign_P1['Bm2'])
ADS_true_sign_P1['Signed_Bp'] = np.log2(ADS_true_sign_P1['Bp1'] / ADS_true_sign_P1['Bp2'])
ADS_true_sign_P2['Signed_Bm'] = np.log2(ADS_true_sign_P2['Bm2'] / ADS_true_sign_P2['Bm1'])
ADS_true_sign_P2['Signed_Bp'] = np.log2(ADS_true_sign_P2['Bp1'] / ADS_true_sign_P2['Bp1'])
ADS_true_signed = pd.concat([ADS_true_sign_P1, ADS_true_sign_P2]).reset_index(drop=True)

# For the minimal model
min_true_sign_P1 = folds_signed_model.copy()
min_true_sign_P1[['Bm1', 'Bm2', 'Bp1', 'Bp2']] = end_true_min[['Bm1', 'Bm2', 'Bp1', 'Bp2']].copy()
min_true_sign_P1['Model'] = 'Minimal'
min_true_sign_P2 = min_sign_P1.copy()
min_true_sign_P1['Signed_Bm'] = np.log2(min_true_sign_P1['Bm1'] / min_true_sign_P1['Bm2'])
min_true_sign_P1['Signed_Bp'] = np.log2(min_true_sign_P1['Bp1'] / min_true_sign_P1['Bp2'])
min_true_sign_P2['Signed_Bm'] = np.log2(min_true_sign_P2['Bm2'] / min_true_sign_P2['Bm1'])
min_true_sign_P2['Signed_Bp'] = np.log2(min_true_sign_P2['Bp1'] / min_true_sign_P2['Bp1'])
min_true_signed = pd.concat([min_true_sign_P1, min_true_sign_P2]).reset_index(drop=True)

# All the signed lo2-fold changes are concatenated, separately for the 'all' and 'true duplicates' sets
signed_all = pd.concat([mixed_signed, nocost_signed, ADS_signed, min_signed]).reset_index(drop=True)
signed_true = pd.concat([mixed_true_signed, nocost_true_signed, ADS_true_signed, min_true_signed]).reset_index(drop=True)

# Finally, the correlations are computed and combined into a dataframe
corr_model = pd.DataFrame(columns=['Model', 'Mut_sigma', 'Mut_ratio', 'Mut_alpha', 'Mut_corr', 'Iter',
                                   'Folder', 'Set', 'Type', 'rho_fold', 'p-val_fold', 'rho_signed',
                                   'p-val_signed', 'rho_ratio', 'p-val_ratio'])
div_correlations = corr_model.copy()

for model in ['Mixed', 'No Cost', 'ADS', 'Minimal']:
    # Correlations when considering all simulated pair
    all_folds_subset = folds_all[folds_all['Model'] == model]
    all_signed_subset = signed_all[signed_all['Model'] == model]
    all_rho_fold = stats.spearmanr(all_folds_subset['Transcription rate'], all_folds_subset['Translation rate'])
    all_rho_signed = stats.spearmanr(all_signed_subset['Signed_Bm'], all_signed_subset['Signed_Bp'])
    all_rho_ratio = stats.spearmanr(all_folds_subset['Protein abundance'], all_folds_subset['Divergence ratio'])

    all_corrs = [model, float(bm_params[1]), ratio, skew_param, correlation, run_name, sim_folder, sim_set, 'All',
                 all_rho_fold[0], all_rho_fold[1], all_rho_signed[0], all_rho_signed[1], all_rho_ratio[0],
                 all_rho_ratio[1]]
    corr_current = corr_model.copy()
    corr_current.loc[0] = all_corrs
    div_correlations = pd.concat([div_correlations, corr_current]).reset_index(drop=True)

    # When considering only the true duplicates
    true_folds_subset = folds_true[folds_true['Model'] == model]
    true_signed_subset = signed_true[signed_true['Model'] == model]
    true_rho_fold = stats.spearmanr(true_folds_subset['Transcription rate'], true_folds_subset['Translation rate'])
    true_rho_signed = stats.spearmanr(true_signed_subset['Signed_Bm'], true_signed_subset['Signed_Bp'])
    true_rho_ratio = stats.spearmanr(true_folds_subset['Protein abundance'], true_folds_subset['Divergence ratio'])

    true_corrs = [model, float(bm_params[1]), ratio, skew_param, correlation, run_name, sim_folder, sim_set,
                  'True duplicates', true_rho_fold[0], true_rho_fold[1], true_rho_signed[0], true_rho_signed[1],
                  true_rho_ratio[0], true_rho_ratio[1]]
    corr_current = corr_model.copy()
    corr_current.loc[0] = true_corrs
    div_correlations = pd.concat([div_correlations, corr_current]).reset_index(drop=True)

# The dataframe containing the correlations is finally exported
div_correlations.reset_index(drop=True).to_csv(os.path.join(path_folder, f'Correlations_{run_name}.csv'), index=False)

# If make_figures has been set to True, they are generated
if make_figs:

    # 1) log2 fold-changes of the three parameters
    if dupli_type == 'all':
        melt_subset = folds_melted[(folds_melted['Model'] != 'WGD couples') & (folds_melted['Model'] != 'SSD couples')]
        hue_order = ['Mixed', 'ADS', 'Minimal', 'No Cost', 'Yeast Data']
        palette = {'Mixed': cm.tab20b.colors[0], 'ADS': cm.tab20b.colors[17], 'Minimal': cm.tab10.colors[9],
                   'Yeast Data': cm.tab20c.colors[4], 'No Cost': cm.tab10.colors[7]}
        box_pairs = [(('Transcription rate', 'Mixed'), ('Translation rate', 'Mixed')),
                     (('Transcription rate', 'Yeast Data'), ('Translation rate', 'Yeast Data')),
                     (('Transcription rate', 'ADS'), ('Transcription rate', 'Yeast Data')),
                     (('Translation rate', 'ADS'), ('Translation rate', 'Yeast Data')),
                     (('Transcription rate', 'Mixed'), ('Transcription rate', 'Yeast Data')),
                     (('Translation rate', 'Mixed'), ('Translation rate', 'Yeast Data')),
                     (('Transcription rate', 'Minimal'), ('Transcription rate', 'Yeast Data')),
                     (('Translation rate', 'Minimal'), ('Translation rate', 'Yeast Data'))]

    elif dupli_type == 'WGD':
        melt_subset = folds_melted[(folds_melted['Model'] != 'Yeast Data') & (folds_melted['Model'] != 'SSD couples')]
        hue_order = ['Mixed', 'ADS', 'Minimal', 'No Cost', 'WGD couples']
        palette = {'Mixed': cm.tab20b.colors[0], 'ADS': cm.tab20b.colors[17], 'Minimal': cm.tab10.colors[9],
                   'WGD couples': cm.tab20c.colors[4], 'No Cost': cm.tab10.colors[7]}
        box_pairs = [(('Transcription rate', 'Mixed'), ('Translation rate', 'Mixed')),
                     (('Transcription rate', 'WGD couples'), ('Translation rate', 'WGD couples')),
                     (('Transcription rate', 'ADS'), ('Transcription rate', 'WGD couples')),
                     (('Translation rate', 'ADS'), ('Translation rate', 'WGD couples')),
                     (('Transcription rate', 'Mixed'), ('Transcription rate', 'WGD couples')),
                     (('Translation rate', 'Mixed'), ('Translation rate', 'WGD couples')),
                     (('Transcription rate', 'Minimal'), ('Transcription rate', 'WGD couples')),
                     (('Translation rate', 'Minimal'), ('Translation rate', 'WGD couples'))]

    elif dupli_type == 'SSD':
        melt_subset = folds_melted[(folds_melted['Model'] != 'Yeast Data') & (folds_melted['Model'] != 'WGD couples')]
        hue_order = ['Mixed', 'ADS', 'Minimal', 'No Cost', 'SSD couples']
        palette = {'Mixed': cm.tab20b.colors[0], 'ADS': cm.tab20b.colors[17], 'Minimal': cm.tab10.colors[9],
                   'SSD couples': cm.tab20c.colors[4], 'No Cost': cm.tab10.colors[7]}
        box_pairs = [(('Transcription rate', 'Mixed'), ('Translation rate', 'Mixed')),
                     (('Transcription rate', 'SSD couples'), ('Translation rate', 'SSD couples')),
                     (('Transcription rate', 'ADS'), ('Transcription rate', 'SSD couples')),
                     (('Translation rate', 'ADS'), ('Translation rate', 'SSD couples')),
                     (('Transcription rate', 'Mixed'), ('Transcription rate', 'SSD couples')),
                     (('Translation rate', 'Mixed'), ('Translation rate', 'SSD couples')),
                     (('Transcription rate', 'Minimal'), ('Transcription rate', 'SSD couples')),
                     (('Translation rate', 'Minimal'), ('Translation rate', 'SSD couples'))]

    changes_only = melt_subset[melt_subset['Parameter'] != 'Divergence ratio']

    plt.figure(figsize=(18, 6))

    violin_folds = sns.violinplot(x='Parameter', y='Log2 fold-change', data=changes_only, hue='Model',
                                  hue_order=hue_order, cut=0, palette=palette)

    statannot.add_stat_annotation(violin_folds, x='Parameter', y='Log2 fold-change', data=changes_only, hue='Model',
                                  hue_order=hue_order, box_pairs=box_pairs,
                                  test='Mann-Whitney', text_format='simple', loc='inside')

    plt.title('Comparison of expression divergence across the four models', size=14)

    folds_fig = violin_folds.get_figure()
    folds_fig.savefig(os.path.join(path_folder, f'folds_violin_{run_name}.png'))
    violin_folds.get_figure().clf()
    print('Violin plots of divergence have been saved.')

    # 2) Divergence ratios across the four models
    ratios_only = melt_subset[melt_subset['Parameter'] == 'Divergence ratio']

    if dupli_type == 'all':
        ratio_pairs = [(('Mixed', 'Mixed'), ('Yeast Data', 'Yeast Data')),
                       (('ADS', 'ADS'), ('Yeast Data', 'Yeast Data')),
                       (('Minimal', 'Minimal'), ('Yeast Data', 'Yeast Data')),
                       (('No Cost', 'No Cost'), ('Yeast Data', 'Yeast Data'))]

    elif dupli_type == 'WGD':
        ratio_pairs = [(('Mixed', 'Mixed'), ('WGD couples', 'WGD couples')),
                       (('ADS', 'ADS'), ('WGD couples', 'WGD couples')),
                       (('Minimal', 'Minimal'), ('WGD couples', 'WGD couples')),
                       (('No Cost', 'No Cost'), ('WGD couples', 'WGD couples'))]

    elif dupli_type == 'SSD':
        ratio_pairs = [(('Mixed', 'Mixed'), ('SSD couples', 'SSD couples')),
                       (('ADS', 'ADS'), ('SSD couples', 'SSD couples')),
                       (('Minimal', 'Minimal'), ('SSD couples', 'SSD couples')),
                       (('No Cost', 'No Cost'), ('SSD couples', 'SSD couples'))]

    box_ratios = sns.boxplot(x='Model', y='Log2 fold-change', order=hue_order, hue_order=hue_order, hue='Model',
                             data=ratios_only, palette=palette)

    statannot.add_stat_annotation(box_ratios, x='Model', y='Log2 fold-change', data=ratios_only, hue='Model',
                                  box_pairs=ratio_pairs, test='Mann-Whitney', text_format='simple', loc='inside')
    box_ratios.get_legend().remove()
    box_ratios.set_ylabel('Log2 divergence ratio')
    box_ratios.set_title('Divergence bias comparison across the four models')

    ratios_fig = box_ratios.get_figure()
    ratios_fig.savefig(os.path.join(path_folder, f'ratios_box_{run_name}.png'))
    box_ratios.get_figure().clf()
    print('Boxplot of divergence ratios has been saved.')

    # 3) Evolution panel
    p_final = pd.read_csv(os.path.join(path_folder, f'Mood_pvalues_{run_name}.csv'))
    p_final.iloc[:, 1:] = np.log10(p_final.iloc[:, 1:])

    fig, axs = plt.subplots(1, 2, figsize=[12, 6])

    sns.lineplot(x='Round', y='Mixed', data=p_final, ax=axs[0], color=cm.tab20b.colors[0], label='Mixed')
    sns.lineplot(x='Round', y='ADS Only', data=p_final, ax=axs[0], color=cm.tab20b.colors[17], label='ADS')
    sns.lineplot(x='Round', y='Minimal', data=p_final, ax=axs[0], color=cm.tab10.colors[9], label='Minimal')
    sns.lineplot(x='Round', y='No Cost', data=p_final, ax=axs[0], color=cm.tab10.colors[7], label='No Cost')

    axs[0].axhline(y=math.log10(0.05), c='k', linestyle='--')
    axs[0].axhline(y=math.log10(0.10), c='k', linestyle='--')

    axs[0].legend()
    axs[0].set_xlabel('Mutation rounds')
    axs[0].set_ylabel("Log10 p-value of Mood's median test")
    axs[0].set_title(f'Simulated protein abundance divergence\ncompared to true divergence ({dupli_type} duplicates)')

    if dupli_loss:
        couples_final = pd.read_csv(os.path.join(path_folder, f'Gene_loss_{run_name}.csv'))
        sns.lineplot(x='Round', y='Mixed', data=couples_final, ax=axs[1], color=cm.tab20b.colors[0], label='Mixed')
        sns.lineplot(x='Round', y='ADS Only', data=couples_final, ax=axs[1], color=cm.tab20b.colors[17], label='ADS')
        sns.lineplot(x='Round', y='Minimal', data=couples_final, ax=axs[1], color=cm.tab10.colors[9], label='Minimal')
        sns.lineplot(x='Round', y='No Cost', data=couples_final, ax=axs[0], color=cm.tab10.colors[7], label='No Cost')
        axs[1].legend()

    axs[1].axhline(y=(0.15 * n_couples), c='k', linestyle='--')

    axs[1].set_xlabel('Mutation rounds')
    axs[1].set_ylabel('Remaining duplicate couples')
    axs[1].set_title('Gene loss through time')

    fig.savefig(os.path.join(path_folder, f'Evolution_panel_{run_name}.png'))
    fig.clf()
    print('Evolution panel has been saved.')

    # 4) Trajectories of each couple through time for the four models
    if dataset == 'Original':
        boundary = 1.1

    elif dataset == 'Corrected':
        boundary = 2.0544

    traj_fig = PdfPages(os.path.join(path_folder, f'Trajectories_{run_name}.pdf'))

    data_mixed = data_all[data_all['Model'] == 'Mixed']
    traj_mixed = evol_funct.trajectories_space(data_mixed, n_couples, total_bm, bm_min, bm_max, bp_min, bp_max,
                                               loss=dupli_loss, model='mixed', boundary=boundary)
    traj_mixed.savefig(traj_fig, format='pdf')
    traj_mixed.clf()

    data_nocost = data_all[data_all['Model'] == 'No Cost']
    traj_nocost = evol_funct.trajectories_space(data_nocost, n_couples, total_bm, bm_min, bm_max, bp_min, bp_max,
                                                loss=dupli_loss, model='no-cost', boundary=boundary)
    traj_nocost.savefig(traj_fig, format='pdf')
    traj_nocost.clf()

    data_ADS = data_all[data_all['Model'] == 'ADS']
    traj_ADS = evol_funct.trajectories_space(data_ADS, n_couples, total_bm, bm_min, bm_max, bp_min, bp_max,
                                             loss=dupli_loss, model='ADS', boundary=boundary)
    traj_ADS.savefig(traj_fig, format='pdf')
    traj_ADS.clf()

    data_min = data_all[data_all['Model'] == 'Minimal']
    traj_min = evol_funct.trajectories_space(data_min, n_couples, total_bm, bm_min, bm_max, bp_min, bp_max,
                                             loss=dupli_loss, model='minimal', boundary=boundary)
    traj_min.savefig(traj_fig, format='pdf')
    traj_min.clf()

    traj_fig.close()

    print('The trajectories panel have been saved.')

    # 5) If full-data saving was enabled, the divergence panels are generated
    if full_data:
        div_panels = PdfPages(os.path.join(path_folder, f'Div_panels_{run_name}.pdf'))

        model_desc = {'Mixed': 'Absolute dosage subfunctionalization with cost-precision tradeoff',
                      'No Cost': 'Absolute dosage subfunctionalization with expression precision constraints',
                      'ADS': 'Absolute dosage subfunctionalization, with ancestral rates defined according to Hausser '
                             'et al. (2019)',
                      'Minimal': 'Absolute dosage subfunctionalization with no relationship between curvature and'
                                 ' ancestral rates'}

        for model in ['Mixed', 'No Cost', 'ADS', 'Minimal']:
            data_panel_current = data_all[data_all['Model'] == model]
            panel_model = evol_funct.divergence_panel(data_panel_current, n_couples, model, model_desc[model],
                                                      loss=dupli_loss)

            panel_model.savefig(div_panels, format='pdf')
            panel_model.clf()
            plt.close()

        div_panels.close()

        print('The divergence panels have been saved.')

    print('All figures have now been saved.')

# End time for the whole script
t_end = time.process_time()

# Generate a file to save the execution time of the script, as well as the parameters of the simulation
time_evol = t1_evol - t0_evol
time_total = t_end - t_start

time_data = [time_total, time_evol, run_name, sim_folder, sim_set, n_couples, bm_params[1], bp_params[1], pop_size,
             skew_param, correlation, bivariate]

time_df = pd.DataFrame(columns=['Time_total', 'Time_evol', 'Iter', 'Folder', 'Set', 'N_pairs', 'Bm_SD', 'Bp_SD',
                                'Population', 'Mut_alpha', 'Mut_corr', 'Bivariate'])

time_df.loc[0] = time_data

time_df.to_csv(os.path.join(path_folder, f'Script_time_{sim_folder}.csv'), mode='a', index=False, header=False)
