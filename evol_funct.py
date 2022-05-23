import scipy.stats as stats
from scipy.stats import (multivariate_normal as mvn,
                         norm)
from scipy.stats._multivariate import _squeeze_output
import pandas as pd
import numpy as np
import scipy.optimize as optimize
import math
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
import statannot

# Useful functions to simulate the divergence of duplicate genes by absolute dosage subfunctionalization under a
# cost-precision tradeoff are defined


# Fold-change between two paralogs
def fold_change(prop_P1, prop_P2, data):
    """Function to calculate log2 fold-change for a property between two duplicates. The two properties are
    provided as dataframe columns, and a dataframe column containing the log2 fold-changes is returned."""

    df = data[[f'{prop_P1}', f'{prop_P2}']].copy()
    df['Fold_Change'] = np.NaN

    for row in range(df.shape[0]):
        value_P1 = df.at[row, f'{prop_P1}']
        value_P2 = df.at[row, f'{prop_P2}']

        if value_P1 == 0 or value_P2 == 0:
            continue

        if value_P1 >= value_P2:
            df.at[row, 'Fold_Change'] = value_P1/value_P2

        elif value_P1 < value_P2:
            df.at[row, 'Fold_Change'] = value_P2/value_P1

    df['Fold_Change'] = np.log2(df['Fold_Change'])

    return df['Fold_Change']


# Divergence ratio between two paralogs when the expression rates are stored as linear values
def div_ratio(bm_P1, bm_P2, bp_P1, bp_P2, data):
    """Function to calculate a linear divergence ratio for a dataframe
    of duplicate couples"""

    df = data[[f'{bm_P1}', f'{bm_P2}', f'{bp_P1}', f'{bp_P2}']].copy()
    df['Bm_ratio'] = np.NaN
    df['Bp_ratio'] = np.NaN
    df['Divergence_ratio'] = np.NaN

    for row in range(df.shape[0]):
        bm1 = df.at[row, f'{bm_P1}']
        bm2 = df.at[row, f'{bm_P2}']

        bp1 = df.at[row, f'{bp_P1}']
        bp2 = df.at[row, f'{bp_P2}']

        if bm1 >= bm2:
            df.at[row, 'Bm_ratio'] = bm1/bm2

        elif bm2 > bm1:
            df.at[row, 'Bm_ratio'] = bm2/bm1

        if bp1 >= bp2:
            df.at[row, 'Bp_ratio'] = bp1/bp2

        elif bp2 > bp1:
            df.at[row, 'Bp_ratio'] = bp2/bp1

    df['Divergence_ratio'] = df['Bm_ratio'] / df['Bp_ratio']

    return df


# Divergence ratio between two paralogs when the expression rates are stored as log10-transformed values
def div_ratio_log(bm_P1, bm_P2, bp_P1, bp_P2, data):
    """Function to calculate a linear divergence ratio for a dataframe
    of duplicate couples"""

    df = data[[f'{bm_P1}', f'{bm_P2}', f'{bp_P1}', f'{bp_P2}']].copy()
    df['Bm_ratio'] = np.NaN
    df['Bp_ratio'] = np.NaN
    df['Divergence_ratio'] = np.NaN

    for row in range(df.shape[0]):
        bm1 = 10**df.at[row, f'{bm_P1}']
        bm2 = 10**df.at[row, f'{bm_P2}']

        bp1 = 10**df.at[row, f'{bp_P1}']
        bp2 = 10**df.at[row, f'{bp_P2}']

        if bm1 >= bm2:
            df.at[row, 'Bm_ratio'] = bm1/bm2

        elif bm2 > bm1:
            df.at[row, 'Bm_ratio'] = bm2/bm1

        if bp1 >= bp2:
            df.at[row, 'Bp_ratio'] = bp1/bp2

        elif bp2 > bp1:
            df.at[row, 'Bp_ratio'] = bp2/bp1

    df['Divergence_ratio'] = df['Bm_ratio'] / df['Bp_ratio']

    return df


# Functions to define bivariate mutational effects distributions. A multivariate skew normal implementation is first
# taken from http://gregorygundersen.com/blog/2020/12/29/multivariate-skew-normal/.
# A mu parameter is added to the constructor, so that the mean of the distribution can be modified (although it will
# probably not be used in the simulations). The possibility to specify a random state is also added
class multivariate_skewnorm:

    def __init__(self, shape, mu, cov=None):
        self.dim = len(shape)
        self.shape = np.asarray(shape)
        self.mean = np.asarray(mu)
        self.cov = np.eye(self.dim) if cov is None else np.asarray(cov)

    def pdf(self, x):
        return np.exp(self.logpdf(x))

    def logpdf(self, x):
        x = mvn._process_quantiles(x, self.dim)
        pdf = mvn(self.mean, self.cov).logpdf(x)
        cdf = norm(0, 1).logcdf(np.dot(x, self.shape))
        return _squeeze_output(np.log(2) + pdf + cdf)

    def rvs_fast(self, size=1, random_state=None):
        aCa = self.shape @ self.cov @ self.shape
        delta = (1 / np.sqrt(1 + aCa)) * self.cov @ self.shape
        cov_star = np.block([[np.ones(1), delta],
                             [delta[:, None], self.cov]])
        x = mvn(np.zeros(self.dim + 1), cov_star).rvs(size=size, random_state=random_state)
        x0, x1 = x[:, 0], x[:, 1:]
        inds = x0 <= 0
        x1[inds] = -1 * x1[inds]
        return x1


# Functions to find the standard deviations of the distribution of mutational effects that respects the desired level
# of mutational bias while also keeping the mean absolute change in protein abundance constant
def diff_prot(dbp, means, bias, correlation, sd_comps, skew, random_state):
    """Computes the mean absolute change of protein abundance for a given combination of standard deviation of
    translational mutational effects and mutational bias for the bivariate case, and compares it to the desired
    magnitude of mean (absolute) mutational effects, set according to a univariate case (independant transcriptional
    and translational mutations). The absolute value of the difference between the two values is returned, so that
    it can be minimized.

    dbp = Standard deviation of the translational mutational effects
    mean = Mean mutational effects. Tuple of form (Mean Bm, Mean Bp).
    bias = Magnitude of mutational bias. A value of 3 means that transcriptional mutations have a mean effect that
           is 3x larger than that of translational mutations.
    correlation = Coefficient of correlation between the transcriptional and translational mutational effects
    sd_comps = Standard deviations of the univariate distributions of mutational effects that the bivariate distribution
               should emulate (provide the same mean absolute change of protein abundance). Tuple of form
               (sd Bm, sd Bp).
    skew = Alpha parameter of the bivariate skew normal distribution (same value for both variables)"""

    # Desired mean absolute change
    bm_dist = stats.skewnorm.rvs(skew, loc=means[0], scale=sd_comps[0], size=100000, random_state=random_state)
    bp_dist = stats.skewnorm.rvs(skew, loc=means[1], scale=sd_comps[1], size=100000, random_state=random_state)
    all_dist = np.concatenate([bm_dist, bp_dist])
    comp = np.mean(np.abs(all_dist)) * 100

    # Obtained mean absolute change
    cov = correlation * math.sqrt((bias * dbp)**2 * (dbp**2))
    cov_mat = np.array([[(bias * dbp)**2, cov], [cov, dbp**2]], dtype=object)

    mut_changes = multivariate_skewnorm((skew, skew), means, cov=cov_mat).rvs_fast(size=100000,
                                                                                   random_state=random_state)

    prot_exp = (1 + mut_changes[:, 0]) * (1 + mut_changes[:, 1])
    prot_percent = np.mean(np.abs(prot_exp - 1)) * 100

    # Difference (to minimize)
    diff = np.abs(prot_percent - comp)

    return diff


def optimize_diff(means, bias, correlation, sd_comps, skew, random_state=None):
    """Function to obtain the standard deviation of translational mutational effects that minimizes the difference with
    the desired mean absolute protein abundance change for a given combination of mutational bias and correlation of
    transcriptional and translational effects. This ensures that varying the level of mutational bias and correlation
    in the bivariate case does not influence the mean protein abundance change of mutations.

    mean = Mean mutational effects. Tuple of form (Mean Bm, Mean Bp)
    bias = Magnitude of mutational bias. A value of 3 means that transcriptional mutations have a mean effect that is
           3x larger than that of translational mutations.
    correlation = Correlation coefficient between transcriptional and translational mutational effects
    sd_comps = Standard deviations of the univariate distributions of mutational effects that the bivariate distribution
               should emulate (provide the same mean absolute change of protein abundance). Tuple of form
               (sd Bm, sd Bp)
    skew = Alpha parameter of the bivariate skew normal distribution (same value for both variables)"""

    opt_dbp = optimize.brute(diff_prot, ((0.01*sd_comps[1], 10*sd_comps[1]),), Ns=10000,
                             args=(means, bias, correlation, sd_comps, skew, random_state))

    dbp = opt_dbp[0]

    return dbp


# Functions to calculate the fitness of specific expression levels
# First, for the ADS + cost-precision mixed model
def variance_dupli(bm1, bm2, p1, p2, alpha_p, cv_0):
    """Function to compute the variance of protein abundance for many duplicate couples at once. The transcription rates
    (bm1 and bm2) and protein abundances (p1 and P2) for each couple are provided as np arrays and an array containing
     the variances for each couple is returned. The variables alpha_p and cv_0 are constants (protein decay rate and
     noise floor) and only need to be provided as floats. It is used in the fit_noise_dupli function below.

    bm1, bm2 = Transcription rate of paralogs P1 and P2, respectively (in mRNA per hour)
    p1, p2 = Protein abundance of paralogs P1 and P2, respectively (in proteins per cell)"""

    with np.errstate(divide='ignore', invalid='ignore'):
        # To avoid warnings raised when p1 == 0 or p2 == 0, since np.where first applies the operation to all elements
        # before testing the condition
        intrinsic_p1 = np.where(p1 == 0, 0, (p1**2) * ((1/p1) + (alpha_p/bm1)))
        intrinsic_p2 = np.where(p2 == 0, 0, (p2**2) * ((1/p2) + (alpha_p/bm2)))

    extrinsic_tot = (np.sqrt((p1**2) * (cv_0**2)) + np.sqrt((p2**2) * (cv_0**2)))**2
    var_tot = intrinsic_p1 + intrinsic_p2 + extrinsic_tot

    return var_tot


def fit_noise_dupli(bm1, bm2, bp1, bp2, pOpt, Q, alpha_m, alpha_p, cv_0):
    """Function to compute the fitness due to protein abundance for a duplicate couple. Protein abundance fluctuations
    du to expression noise are taken into account, with intrinsic noise being completely decoupled between the two
    paralogs, while extrinsic noise is perfectly correlated between them. The transcription rates of both genes (bm1
    and bm2), their protein abundance (p1 and p2), the optimal cumulative protein abundance (pOpt) and the noise
    sensitivity (Q, curvature of the fitness function) are used. The function is assumed to be parabolic, with an
    optimum at (pOpt, 0.42), where fitness is in 1/hour.

    bm1, bm2 = Transcription rate of paralogs P1 and P2, respectively (in mRNAs per hour)
    bp1, bp2 = Translation rate of paralogs P1 and P2, respectively (in proteins per mRNA per hour)
    pOpt = Optimum of the fitness function of cumulative protein abundance (in proteins per cell)
    Q = Noise sensitivity, as defined by Hausser et al., 2019
    alpha_m = mRNA decay rate (constant, in 1/hour)
    alpha_p = protein decay rate (constant, in 1/hour)
    cv_0 = Noise floor (minimal coefficient of variation of protein abundance), as reported by Hausser et al., 2019"""

    f_d2 = -2*(Q/pOpt)

    p1 = (bm1 * bp1) / (alpha_m * alpha_p)
    p2 = (bm2 * bp2) / (alpha_m * alpha_p)

    # Variance of the protein abundance distribution
    var_tot = variance_dupli(bm1, bm2, p1, p2, alpha_p, cv_0)

    # Parameters of the parabola of fitness according to protein abundance
    a = f_d2/2
    b = -2 * a * pOpt
    c = 0.42 - (a * pOpt ** 2) - (b * pOpt)

    # Mean fitness when expression noise is taken into account
    # E(Fit) = a * E(x**2) + b * E(x) + c
    mean_x2 = var_tot + (p1 + p2)**2
    mean_fit = a*mean_x2 + b*(p1 + p2) + c

    return mean_fit


def fit_global_dupli(bm1, bm2, bp1, bp2, pOpt, Q, alpha_m, alpha_p, cv_0, lm, c_m):
    """Function to compute the fitness (in 1/hour) of a duplicate couple for a given protein abundance,
    taking into account expression noise and transcription costs.

    bm1, bm2 = Transcription rate of paralogs P1 and P2, respectively (in mRNAs per hour)
    bp1, bp2 = Translation rate of paralogs P1 and P2, respectively (in proteins per mRNA per hour
    pOpt = Optimal cumulative protein abundance (in proteins per cell)
    Q = Noise sensitivity, as defined by Hausser et al., 2019
    alpha_m = mRNA decay rate (constant, in 1/hour)
    alpha_p = protein decay rate (constant, in 1/hour)
    cv_0 = Noise floor (minimal coefficient of variation of protein abundance), as reported by Hausser et al., 2019
    lm = Length of the mRNA, assumed to remain equal for both duplicates (in nucleotides)
    c_m = Fitness cost of transcription per nucleotide, as calculated by Hausser et al., 2019 (in 1/hours)"""

    fit_abun = fit_noise_dupli(bm1, bm2, bp1, bp2, pOpt, Q, alpha_m, alpha_p, cv_0)

    exp_cost = c_m * lm * (bm1 + bm2)

    fit_glob = fit_abun - exp_cost

    return fit_glob


def singleton_opt(params, pOpt, Q, lm, c_m, alpha_p, cv_0):
    """Fitness function that will be used to optimize (transcription rate and mean protein abundance)
    the ancestral singleton

    params = Tuple (current protein abundance, current transcription rate) in units (proteins per cell, MRNAs per hour)
    pOpt = Optimal protein abundance (in proteins per cell)
    Q = Noise sensitivity, as defined by Hausser et al., 2019
    lm = mRNA length (in nucleotides)
    c_m = Fitness cost of transcription per nucleotide, as calculated by Hausser et al., 2019 (in 1/hours)
    alpha_p = Protein decay rate (in 1/hours)
    cv_0 = Noise floor (minimal possible coefficient of variation for protein abundance),
     as reported by Hausser et al., 2019"""

    x, bm = params

    # Variance of the protein abundance distribution
    var = (x**2) * ((1/x) + (alpha_p/bm) + (cv_0**2))

    # Parameters of the parabola of fitness
    f_d2 = -2 * (Q / pOpt)
    a = f_d2/2
    b = -2 * a * pOpt
    c = 0.42 - (a * pOpt ** 2) - (b * pOpt)

    # Mean fitness due to protein abundance and global fitness taking transcription cost into account
    mean_x2 = var + x**2
    mean_fit = a*mean_x2 + b*x + c

    exp_cost = c_m * lm * bm
    fit_glob = -1*(mean_fit - exp_cost)

    return fit_glob


def ancestral_opt(params, pOpt, Q, lm, c_m, alpha_m, alpha_p, cv_0):
    """Fitness function that will be used to optimize (transcription and translation rates)
    the ancestral singleton

    params = Tuple (current transcription rate, current translation rate) in units (mRNAs per hour, proteins per mRNA
    per hour)
    pOpt = Optimal protein abundance (in proteins per cell)
    Q = Noise sensitivity, as defined by Hausser et al., 2019
    lm = mRNA length (in nucleotides)
    c_m = Fitness cost of transcription per nucleotide, as calculated by Hausser et al., 2019 (in 1/hours)
    alpha_m = mRNA decay rate (in 1/hours)
    alpha_p = Protein decay rate (in 1/hours)
    cv_0 = Noise floor (minimal possible coefficient of variation for protein abundance),
    as reported by Hausser et al., 2019"""

    bm, bp = params

    x = (bm * bp) / (alpha_m * alpha_p)

    # Variance of the protein abundance distribution

    var = np.where(x == 0, 0, ((x**2) * ((1/x) + (alpha_p/bm) + (cv_0**2))))

    # Parameters of the parabola of fitness
    f_d2 = -2 * (Q / pOpt)
    a = f_d2/2
    b = -2 * a * pOpt
    c = 0.42 - (a * pOpt ** 2) - (b * pOpt)

    # Mean fitness due to protein abundance and global fitness taking transcription cost into account
    mean_x2 = var + x**2
    mean_fit = a*mean_x2 + b*x + c

    exp_cost = c_m * lm * bm
    fit_glob = -1*(mean_fit - exp_cost)

    return fit_glob


# Still for the mixed model, a variant is also defined where protein abundance fluctuations are completely independent
# between the two duplicates (noise is not divided into intrinsic and extrinsic components)
def variance_ind(bm1, bm2, p1, p2, alpha_p, cv_0):
    """Total variance of protein abundance for a duplicate couple, if protein abondance fluctuations are considered
    totally independent between the two copies. Designed to receive bm1, bm1, p1 and p2 as numpy arrays and to return
    an array of variances.

    bm1, bm2 = Transcription rate of paralogs P1 and P2, respectively (in mRNAs per hour)
    p1, p2 = Protein abundance for paralogs P1 and P2, respectively (in proteins per cell)
    alpha_p = Protein decay rate (in 1/hour)
    cv_0 = Noise floor (minimal coefficient of variation of protein abundance), as reported by Hausser et al., 2019"""

    with np.errstate(divide='ignore', invalid='ignore'):
        # To avoid the warnings raised when bm1 == 0 or bm2 == 0, since np.where first applies the operation to all
        # elements before testing the condition
        var_p1 = np.where(bm1 == 0, 0, p1**2 * ((1/p1) + (alpha_p/bm1) + cv_0**2))
        var_p2 = np.where(bm2 == 0, 0, p2**2 * ((1/p2) + (alpha_p/bm2) + cv_0**2))

    var_tot = var_p1 + var_p2

    return var_tot


def fit_noise_ind(bm1, bm2, bp1, bp2, pOpt, Q, alpha_m, alpha_p, cv_0):
    """Function to compute the fitness due to protein abundance while taking expression noise into account for a
    duplicate couple when protein abundance fluctuations are assumed to be totally independent between the two copies.
    The fitness function is assumed to be a parabola of vertex (pOpt, 0.42), 1/hours, with a curvature according the
    noise sensitivity Q. The bm and bp rates are given as numpy arrays, so that the function returns one array of
    corresponding fitness values.

    bm1, bm2 = Transcription rate of paralogs P1 and P2, respectively (in mRNAs per hour)
    bp1, bp2 = Translation rate of paralogs P1 and P2, respectively (in proteins per mRNA per hour)
    pOpt = Optimal cumulative protein abundance (in proteins per cell)
    Q = Noise sensitivity, as defined by Hausser et al., 2019
    alpha_m = mRNA decay rates (in 1/hour)
    alpha_p = Protein decay rates (in 1/hour)
    cv_0 = Noise floor (minimal coefficient of variation of protein abundance), as reported by Hausser et al., 2019."""

    f_d2 = -2 * (Q / pOpt)

    p1 = (bm1 * bp1) / (alpha_m * alpha_p)
    p2 = (bm2 * bp2) / (alpha_m * alpha_p)

    # Variance of the protein abundance distribution
    var_tot = variance_ind(bm1, bm2, p1, p2, alpha_p, cv_0)

    # Parameters of the parabola of fitness according to protein abundance
    a = f_d2 / 2
    b = -2 * a * pOpt
    c = 0.42 - (a * pOpt ** 2) - (b * pOpt)

    # Mean fitness when expression noise is taken into account
    # E(Fit) = a * E(x**2) + b * E(x) + c
    mean_x2 = var_tot + (p1 + p2) ** 2
    mean_fit = a * mean_x2 + b * (p1 + p2) + c

    return mean_fit


def fit_global_ind(bm1, bm2, bp1, bp2, pOpt, Q, alpha_m, alpha_p, cv_0, lm, c_m):
    """Function to compute the fitness (in 1/hour) of a duplicate couple for a given protein abundance,
    taking into account expression noise and transcription costs, under the assumption that protein abundance
    fluctutations are totally independent between the two copies.

    bm1, bm2 = Transcription rate of paralogs P1 and P2, respectively (in mRNAs per hour)
    bp1, bp2 = Translation rate of paralogs P1 and P2, respectively (in proteins per mRNA per hour
    pOpt = Optimal cumulative protein abundance (in proteins per cell)
    Q = Noise sensitivity, as defined by Hausser et al., 2019
    alpha_m = mRNA decay rate (constant, in 1/hour)
    alpha_p = protein decay rate (constant, in 1/hour)
    cv_0 = Noise floor (minimal coefficient of variation of protein abundance), as reported by Hausser et al., 2019
    lm = Length of the mRNA, assumed to remain equal for both duplicates (in nucleotides)
    c_m = Fitness cost of transcription per nucleotide, as calculated by Hausser et al., 2019 (in 1/hours)"""

    fit_abun = fit_noise_ind(bm1, bm2, bp1, bp2, pOpt, Q, alpha_m, alpha_p, cv_0)

    exp_cost = c_m * lm * (bm1 + bm2)

    fit_glob = fit_abun - exp_cost

    return fit_glob


# Second, strictly for absolute dosage subfunctionalization (ADS)
def fit_parabola(bm1, bm2, bp1, bp2, pOpt, Q, alpha_m, alpha_p):
    """Function returning fitness (in 1/hour) for a duplicate couple depending solely on mean cumulative protein
    abundance. The fitness function is assumed to be a parabola of optimum (pOpt, 0.42), with a curvature set by noise
    sensitivity Q.

    bm1, bm2 = Transcription rate of paralogs P1 and P2, respectively (in mRNAs per hour)
    bp1, bp2 = Translation rate of paralogs P1 and P2, respectively (in proteins per mRNA)
    pOpt = Optimal cumulative protein abundance (in proteins per cell)
    Q = Noise sensitivity, as described by Hausser et al., 2019
    alpha_m = mRNA decay rate (in 1/hour)
    alpha_p = Protein decay rate (in 1/hour)"""

    # Protein abundances
    p1 = (bm1 * bp1) / (alpha_m * alpha_p)
    p2 = (bm2 * bp2) / (alpha_m * alpha_p)
    exp_tot = p1 + p2

    # Parabola of fitness
    f_d2 = -2 * (Q / pOpt)
    a = f_d2/2
    b = -2 * a * pOpt
    c = 0.42 - (a * pOpt ** 2) - (b * pOpt)

    fit = a * exp_tot**2 + b * exp_tot + c

    return fit


def fit_para_single(bm, bp, pOpt, Q, alpha_m, alpha_p):
    """Function returning fitness (in 1/hour) for a singleton gene according to protein abundance
    (without taking expression noise or transcription costs into account). The function is assumed
    to be a parabola with vertex (pOpt, 0.42) and a curvature set by the noise sensitivity Q.

    bm = Transcription rate of the singleton gene (in mRNAs per hour)
    bp = Translation rate of the singleton gene (in proteins per mRNA per hour)
    pOpt = Optimal protein abundance (in proteins per cell)
    Q = Noise sensitivity, as defined by Hausser et al., 2019
    alpha_m = mRNA decay rate (in 1/hour)
    alpha_p = Protein decay rate (in 1/hour)"""

    # Protein abundance
    p = (bm * bp) / (alpha_m * alpha_p)

    # Parabola of fitness
    f_d2 = -2 * (Q / pOpt)
    a = f_d2/2
    b = -2 * a * pOpt
    c = 0.42 - (a * pOpt ** 2) - (b * pOpt)

    fit = a * p**2 + b * p + c

    return fit


# Third, for the cost-precision model (cost-precision tradeoff without any absolute dosage subfunctionalization)
def fit_noise_single(bm, bp, pOpt, Q, alpha_m, alpha_p, cv_0):
    """Function to compute the fitness (in 1/hour) from protein abundance for a singleton gene, taking into account
    expression noise and transcription costs. The fitness function is assumed to be a parabola of vertex (pOpt, 0.42),
    with a curvature set by noise sensitivity Q.

     bm = Transcription rate of the gene (in mRNAs per hour)
     bp = Translation rate of the gene (in proteins per mRNA per hour)
     pOpt = Optimal protein abundance (in proteins per cell)
     Q = Noise sensitivity, as defined by Hausser et al., 2019
     alpha_m = mRNA decay rate (in 1/hour)
     alpha_p = protein decay rate (in 1/hour)
     cv_0 = Noise floor (minimal coefficient of variation of protein abundance), as reported by Hausser et al., 2019"""

    p = (bm * bp) / (alpha_m * alpha_p)

    # Variance of the protein abundance distribution
    var = np.where(p == 0, 0, (p**2) * ((1/p) + (alpha_p/bm) + (cv_0**2)))

    # Parameters of the parabola of fitness
    f_d2 = -2 * (Q / pOpt)
    a = f_d2/2
    b = -2 * a * pOpt
    c = 0.42 - (a * pOpt ** 2) - (b * pOpt)

    # Mean fitness due to protein abundance
    mean_x2 = var + p**2
    mean_fit = a*mean_x2 + b*p + c

    return mean_fit


def fit_glob_single(bm, bp, pOpt, Q, alpha_m, alpha_p, cv_0, lm, c_m):
    """Function to compute the fitness of a singleton gene according to the cost-precision trade-off.

    bm = Transcription rate (in mRNAs per hour)
    bp = Translation rate (in proteins per mRNA per hour)
    pOpt = Optimal protein abundance (in proteins per cell)
    Q = Noise sensitivity for the fitness function, as defined by Hausser et al., 2019
    alpha_m = mRNA decay rate (in 1/hour)
    alpha_p = Protein decay rate (in 1/hour)
    cv_0 = Noise floor (minimal coefficient of variation of protein abundance), as reported by Hausser et al., 2019
    lm = Length of the mRNA (in nucleotides)
    c_m = Fitness cost of transcription per nucleotide, as calculated by Hausser et al., 2019 (in 1/hours)"""

    fit_exp = fit_noise_single(bm, bp, pOpt, Q, alpha_m, alpha_p, cv_0)

    exp_cost = lm * c_m * bm

    fit_glob = fit_exp - exp_cost

    return fit_glob


def fit_glob_2ind(bm1, bm2, bp1, bp2, pOpt, Q, alpha_m, alpha_p, cv_0, lm, c_m):
    """Function to compute the fitness according to the expression level of two duplicates when selection acts on
    their individual expression levels rather than on cumulative protein abundance. Both paralogs have the same
    optimal mean protein abundance (pOpt), which is set according to that of the ancestral singleton as well as the
    duplication effect. Expression noise as well as transcription costs are taken into account.

    bm1, bm2 = Transcription rate of paralogs P1 and P2, respectively (in mRNAs per hour)
    bp1, bp2 = Translation rate of paralogs P1 and P2, respextively (in proteins per mRNA per hour)
    pOpt = Optimal protein abundance for each duplicates, assumed to be equal (in proteins per cell)
    Q = Noise sensitivity for the two fitness functions (assumed equal), as defined by Hausser et al., 2019
    alpha_m = mRNA decay rate (in 1/hour)
    alpha_p = Protein decay rate (in 1/hour)
    cv_0 = Noise floor (minimal coefficient of variation of protein abundance), as reported by Hausser et al., 2019
    lm = Length of the mRNA, assumed to remain equal for both duplicates (in nucleotides)
    c_m = Fitness cost of transcription per nucleotide, as calculated by Hausser et al., 2019 (in 1/hours)"""

    fit_exp1 = fit_noise_single(bm1, bp1, pOpt, Q, alpha_m, alpha_p, cv_0)/0.42
    fit_exp2 = fit_noise_single(bm2, bp2, pOpt, Q, alpha_m, alpha_p, cv_0)/0.42

    fit_mean = (fit_exp1 * fit_exp2) * 0.42

    exp_cost = lm * c_m * (bm1 + bm2)

    fit_glob = fit_mean - exp_cost

    return fit_glob


# Functions to compute the fixation probability of mutations
def metropolis(fi, fj, N):
    """Function to compute fixation probability of a mutation according to the Metropolis criterion,
     where fj is the new fitness and fi is the former fitness.

     fi = Ancestral fitness (before the mutation)
     fj = Mutant fitness
     N = Population size"""

    # In case there are negative values (which should not happen)
    fi = np.where(fi <= 0, 1e-60, fi)
    fj = np.where(fj <= 0, 1e-60, fj)

    fi = fi.astype('float64')
    fj = fj.astype('float64')

    with np.errstate(all='raise'):
        delta_fit = np.log10(fi) - np.log10(fj)

    with np.errstate(over='ignore'):
        # To avoid the warnings due to overflows when delta_fit < 0 (because np.where first applies the operation to
        # all elements before testing the condition
        prob = np.where(delta_fit < 0, 1, np.exp(-2*N*delta_fit))

    return prob


def sella_hirsh(fi, fj, Ne):
    """Function to compute fixation probability of a mutation according to a Moran birth-death process. From
    Sella & Hirsh, 2005. The former fitness is fi and the new fitness is fj.

    fi = Ancestral fitness (before the mutation)
    fj = Mutant fitness
    Ne = Effective population size"""

    fit_ratio = fi/fj

    prob = np.where(fit_ratio == 1, 1/Ne, (1 - (fi/fj)) / (1 - (fi/fj)**Ne))

    return prob


def kimura(fi, fj, Ne):
    """Function to compute the fixation probability of a mutation according to Kimura's diffusion approximation.
    The former fitness fi and the mutant fitness fj are given as numpy arrays. An haploid population is assumed.

    fi = Ancestral fitness (before the mutation)
    fj = Mutant fitness
    Ne = Effective population size"""

    s = (fj / fi) - 1

    prob = np.where(s == 0, 1/Ne, (1 - np.exp(-s)) / (1 - np.exp(-2 * Ne * s)))

    return prob


# Functions to use during the evolution of the duplicates
def optimize_single(pOpt_anc, Q, alpha_m, alpha_p, cv_0, c_m, lm):
    """Returns optimal bm and bp rates for an ancestral singleton according to a cost-precision tradeoff. As the
    magnitude of noise depends on the expression level, a mean protein abundance slightly lower than the optimal protein
    abundance is optimal, so the rates can not simply be calculated from the formula given by Hausser et al. to obtain
    the optimal bp/bm ratio. The optimization is performed using a grid-based 'brute force' approach. Both values
    are given as floats, and two floats are returned (does not support numpy arrays)

    pOpt_anc = Optimal protein abundance (in proteins per cell)
    Q = Noise sensitivity, as defined by Hausser et al., 2019
    alpha_m = mRNA decay rate (in 1/hour)
    alpha_p = Protein decay rate (in 1/hour)
    cv_0 = Noise floor (minimal coefficient of variation of protein abundance), as reported by Hausser et al., 2019
    c_m = Fitness cost of transcription per nucleotide, as calculated by Hausser et al., 2019 (in 1/hours)
    lm = mRNA length (in nucleotides)"""

    pOpt_bounds = ((pOpt_anc - (0.5 * pOpt_anc)), (pOpt_anc + (0.5 * pOpt_anc)))

    # Optimal transcription rate according to Hausser et al. as a starting point
    bm_est = math.sqrt(pOpt_anc * alpha_p * alpha_m * (Q / (c_m * lm * alpha_m)))
    bm_bounds = ((bm_est - (0.5 * bm_est)), (bm_est + (0.5 * bm_est)))

    bm_res = []
    prot_res = []

    for i in range(5):
        opt_single = optimize.brute(singleton_opt, (pOpt_bounds, bm_bounds), args=(pOpt_anc, Q, lm, c_m, alpha_p,
                                                                                   cv_0), Ns=2000)
        bm_res.append(opt_single[1])
        prot_res.append(opt_single[0])

    bm_opti = np.mean(bm_res)
    prot_opti = np.mean(prot_res)

    bp_opti = (prot_opti * alpha_m * alpha_p) / bm_opti

    return bm_opti, bp_opti


def gene_loss(ancestral_rates, ancestral_fit, fit_funct, pop, args_fit=()):
    """Function to perform gene loss if it is neutral (dFit >= - 1/Ne) at any point in the simulation. A numpy array
    of the ancestral states (prior to the current mutation) of all simulated duplicate couples is given, and a new
    version of this array where all allowed loss-of-function mutations have been done is returned. The number of
    intact duplicate couples is also returned

    ancestral_rates = Numpy array containing the transcription and translation rates for paralogs P1 and P2 of each
                      simulated duplicate couple. Each row of the array is a distinct couple, while its columns are,
                      in order: [Couple (int identifier), bm_P1, bp_P1, bm_P2, bp_P2]
    ancestral_fit = Numpy array containing the fitness values calculated for each of the simulated couples at the
                   start of the current mutation round (ancestral fitness). This array is unidimensional
    fit_funct = Fitness function to be called when computing the fitness effect of gene loss. Its first four arguments
                need to be bm_P1, bm_P2, bp_P1 and bp_P2, in this order.
    pop = Population size assumed in the current run of the simulation.
    args_fit = Tuple of all arguments required by fit_funct, other than bm_P1, bm_P2, bp_P1 and bp_P2"""

    loss_array = ancestral_rates.copy()  # Copy of the array that will be modified and returned
    zeros_val = np.zeros(ancestral_rates.shape[0])

    # Fitness effects of gene loss, both for P1 and P2
    diff_P1 = ancestral_fit - fit_funct(zeros_val, ancestral_rates[:, 3], zeros_val, ancestral_rates[:, 4], *args_fit)
    diff_P2 = ancestral_fit - fit_funct(ancestral_rates[:, 1], zeros_val, ancestral_rates[:, 2], zeros_val, *args_fit)

    # Testing whether gene loss is neutral, according to the fitness effects computed above
    diff_P1 = np.where(diff_P1 < (1 / pop), 0, 1)
    diff_P2 = np.where(diff_P2 < (1 / pop), 0, 1)

    # Identifying the least expressed copy of each couple
    # False, thus 0, when the copy of interest (left) is the least expressed
    P1_low = (ancestral_rates[:, 1] * ancestral_rates[:, 2]) >= (ancestral_rates[:, 3] * ancestral_rates[:, 4])
    P2_low = (ancestral_rates[:, 3] * ancestral_rates[:, 4]) >= (ancestral_rates[:, 1] * ancestral_rates[:, 2])

    # Identifying copies that are the least expressed while their loss is neutral
    # Only the positions where both conditions are True will remain as 0
    P1_del = diff_P1 + P1_low
    P2_del = diff_P2 + P2_low

    # Gene loss is performed when appropriate
    loss_array[:, 1:3] = loss_array[:, 1:3] * (np.where(P1_del > 0, 1, 0)[:, np.newaxis])
    loss_array[:, 3:5] = loss_array[:, 3:5] * (np.where(P2_del > 0, 1, 0)[:, np.newaxis])

    # The number of intact duplicate couples is measured
    n_dupli = ancestral_rates.shape[0] - (np.where(loss_array[:, 1] == 0, 1, 0).sum() +
                                          np.where(loss_array[:, 3] == 0, 1, 0).sum())

    return loss_array, n_dupli


def muta_bidirectional(current_rates, n_couples, bm_limits, bp_limits, bm_to_low, bm_mid, bm_to_high, bp_to_low, bp_mid,
                       bp_to_high, P1_choice):
    """Function to produce the final array of mutations when genes have different mutational effects distributions
    (varying directions of asymmetry) depending on their current rates of transcription and translation.

    current_rates = Array of current transcription and translation rates (shape = (n_couples, 5), columns = [Couple,
                    Bm1, Bp1, Bm2, Bp2])
    n_couples = Number of duplicate couples included in the current simulation
    bm_limits = Tuple of (min, max) values of transcription rates past which asymmetry is introduced
    bp_limits = Tuple of (min, max) values of translation rates past which asymmetry is introduced
    All the other arguments are arrays of shape (n_couples, )
    bm_to_low = Distribution of transcriptional mutational effects biases towards a reduction
    bm_mid = Distribution of transcriptional mutational effects with no bias
    bm_to_high = Distribution of transcriptional mutational effects biased towards an increase
    bp_to_low = Distribution of translational mutational effects biased towards a reduction
    bp_mid = Distribution of translational mutational effects with no bias
    bp_to_high = Distribution of translational mutational effects biased towards an increase
    P1_choice = Array of shape (n_couples, ) identifying where paralog P1 will receive a mutation (1 for a mutation,
                0 for no mutation)"""

    # Making the transcriptional component of the array of mutations
    low_bm_P1 = np.where(current_rates[:, 1] <= bm_limits[0], 1, 0)
    low_bm_P2 = np.where(current_rates[:, 3] <= bm_limits[0], 1, 0)
    high_bm_P1 = np.where(current_rates[:, 1] >= bm_limits[1], 1, 0)
    high_bm_P2 = np.where(current_rates[:, 3] >= bm_limits[1], 1, 0)
    mid_bm_P1 = np.where((current_rates[:, 1] > bm_limits[0]) & (current_rates[:, 1] < bm_limits[1]), 1, 0)
    mid_bm_P2 = np.where((current_rates[:, 3] > bm_limits[0]) & (current_rates[:, 3] < bm_limits[1]), 1, 0)

    # Lower values will be biased towards an increase, while higher values will be biased towards a decrease
    bm_to_low_P1 = high_bm_P1 * bm_to_low
    bm_to_low_P2 = high_bm_P2 * bm_to_low
    bm_to_high_P1 = low_bm_P1 * bm_to_high
    bm_to_high_P2 = low_bm_P2 * bm_to_high
    bm_mid_P1 = mid_bm_P1 * bm_mid
    bm_mid_P2 = mid_bm_P2 * bm_mid

    bm_P1 = bm_to_low_P1 + bm_to_high_P1 + bm_mid_P1
    bm_P2 = bm_to_low_P2 + bm_to_high_P2 + bm_mid_P2

    # Similarly producing the translational component of the final array of mutations
    low_bp_P1 = np.where(current_rates[:, 2] <= bp_limits[0], 1, 0)
    low_bp_P2 = np.where(current_rates[:, 4] <= bp_limits[0], 1, 0)
    high_bp_P1 = np.where(current_rates[:, 2] >= bp_limits[1], 1, 0)
    high_bp_P2 = np.where(current_rates[:, 4] >= bp_limits[1], 1, 0)
    mid_bp_P1 = np.where((current_rates[:, 2] > bp_limits[0]) & (current_rates[:, 2] < bp_limits[1]), 1, 0)
    mid_bp_P2 = np.where((current_rates[:, 4] > bp_limits[0]) & (current_rates[:, 4] < bp_limits[1]), 1, 0)

    bp_to_low_P1 = high_bp_P1 * bp_to_low
    bp_to_low_P2 = high_bp_P2 * bp_to_low
    bp_to_high_P1 = low_bp_P1 * bp_to_high
    bp_to_high_P2 = low_bp_P2 * bp_to_high
    bp_mid_P1 = mid_bp_P1 * bp_mid
    bp_mid_P2 = mid_bp_P2 * bp_mid

    bp_P1 = bp_to_low_P1 + bp_to_high_P1 + bp_mid_P1
    bp_P2 = bp_to_low_P2 + bp_to_high_P2 + bp_mid_P2

    # Decision on whether P1 or P2 is mutated
    P2_choice = np.abs(P1_choice - 1)

    # Final array of mutational effects
    mutations = np.zeros((n_couples, 4))
    mutations[:, 0] = bm_P1 * P1_choice
    mutations[:, 1] = bp_P1 * P1_choice
    mutations[:, 2] = bm_P2 * P2_choice
    mutations[:, 3] = bp_P2 * P2_choice

    return mutations


# Plotting functions and necessary data cleaning functions
def panel_ancestral(bm_mixed, bp_mixed, bm_ADS, bp_ADS, pOpt_mixed, pOpt_ADS, Q_mixed, Q_ADS, alpha_m,
                    alpha_p, cv_0, lm, c_m):
    """Function to produce a 4-panel figure validating the optimality of two different ancestral singletons,
     respectively used for simulations according to the mixed model (ADS + cost-precision) and to the ADS only model.
     Plots of fitness across transcription rate and translation rate variations are made for both cases.

    bm_mixed, bp_mixed = Optimal (ancestral) transcription and translation rates according to the mixed model
                         (respectively in mRNAs per hour and in proteins per mRNA per hour
    bm_ADS, bp_ADS = Optimal (ancestral) transcription and translation rates according to ADS only (respectively in
                     mRNAs per hour and in proteins per mRNa per hour
    pOpt_mixed, pOpt_ADS = Optimal (ancestral) protein abundance, which will be equal for the two singletons
                           considered (in proteins per cell)
    Q_mixed, Q_ADS = Noise sensitivity, as defined by Hausser et al., 2019
    alpha_m = mRNA decay rate (in 1/hour)
    alpha_p = Protein decay rate (in 1/hour)
    cv_0 = Noise floor (minimal coefficient of variation of protein abundance), as reported by Hausser et al., 2019
    lm = pre-mRNA length, only used for the calculation of fitness according to the mixed model
    c_m = Expression cost per transcribed nucleotide (Hausser et al., 2019), only used for calculations according to
          the mixed model"""

    dbm_mixed = np.linspace(-0.95 * bm_mixed, 0.95 * bm_mixed, 250)
    dbp_mixed = np.linspace(-0.95 * bp_mixed, 0.95 * bp_mixed, 250)
    dbm_ADS = np.linspace(-0.95 * bm_ADS, 0.95 * bm_ADS, 250)
    dbp_ADS = np.linspace(-0.95 * bp_ADS, 0.95 * bp_ADS, 250)

    bms_mixed = dbm_mixed + bm_mixed
    bps_mixed = dbp_mixed + bp_mixed
    bms_ADS = dbm_ADS + bm_ADS
    bps_ADS = dbp_ADS + bp_ADS

    prot_mixed_bm = (bms_mixed * bp_mixed) / (alpha_m * alpha_p)
    prot_mixed_bp = (bm_mixed * bps_mixed) / (alpha_m * alpha_p)

    fit_bm_mixed = -1 * singleton_opt((prot_mixed_bm, bms_mixed), pOpt_mixed, Q_mixed, lm, c_m, alpha_p, cv_0) - 0.42
    fit_bp_mixed = -1 * singleton_opt((prot_mixed_bp, bm_mixed), pOpt_mixed, Q_mixed, lm, c_m, alpha_p, cv_0) - 0.42
    fit_bm_ADS = fit_para_single(bms_ADS, bp_ADS, pOpt_ADS, Q_ADS, alpha_m, alpha_p) - 0.42
    fit_bp_ADS = fit_para_single(bm_ADS, bps_ADS, pOpt_ADS, Q_ADS, alpha_m, alpha_p) - 0.42

    fig, axs = plt.subplots(2, 2, figsize=[16, 16])
    plasma = cm.get_cmap('plasma', 4)  # Colors for the lines

    axs[0, 0].plot(dbm_mixed, fit_bm_mixed, c=plasma.colors[0])
    axs[0, 0].axvline(0, 0, 0.95, c='k', linestyle='--')
    axs[0, 0].set_title('Transcription rate variation, mixed model')
    axs[0, 0].set_ylabel('Fitness difference')
    axs[0, 0].set_xlabel('Transcription rate change (mRNA per hour)')

    axs[0, 1].plot(dbp_mixed, fit_bp_mixed, c=plasma.colors[1])
    axs[0, 1].axvline(0, 0, 0.95, c='k', linestyle='--')
    axs[0, 1].set_title('Translation rate variation, mixed model')
    axs[0, 1].set_ylabel('Fitness difference')
    axs[0, 1].set_xlabel('Translation rate change (protein per mRNA per hour)')

    axs[1, 0].plot(dbm_ADS, fit_bm_ADS, c=plasma.colors[2])
    axs[1, 0].axvline(0, 0, 0.95, c='k', linestyle='--')
    axs[1, 0].set_title('Transcription rate variation, ADS model')
    axs[1, 0].set_ylabel('Fitness difference')
    axs[1, 0].set_xlabel('Transcription rate change (mRNA per hour)')

    axs[1, 1].plot(dbp_ADS, fit_bp_ADS, c=plasma.colors[3])
    axs[1, 1].axvline(0, 0, 0.95, c='k', linestyle='--')
    axs[1, 1].set_title('Translation rate variation, ADS model')
    axs[1, 1].set_ylabel('Fitness difference')
    axs[1, 1].set_xlabel('Translation rate change (protein per mRNA per hour)')

    return fig


def process_data(sim_data, loss=False):
    """Function that processes the simulation data into separate dataframe for duplicate couples and for singletons

    sim_data = 'data' dataframe saved at the end of the simulation
    loss = Boolean indicating whether gene loss was allowed during the simulation run used"""

    div_data = sim_data.copy()
    div_data.iloc[:, 2:] = div_data.iloc[:, 2:].where(div_data.iloc[:, 2:] != np.inf, other=0)

    if loss:
        single_df = div_data[(div_data['Bm1'] == 0) | (div_data['Bm2'] == 0)].copy()
        duplicates = div_data[(div_data['Bm1'] != 0) & (div_data['Bm2'] != 0)].copy()

        # The singletons df is cleaned up (to remove zeros)
        single_P1 = single_df.iloc[:, np.r_[0:4, 6, 8]].copy()
        single_P2 = single_df.iloc[:, np.r_[0:2, 4, 5, 7, 9]].copy()

        single_P1.columns = ['Round', 'Couple', 'Bm', 'Bp', 'Prot', 'cv']
        single_P2.columns = ['Round', 'Couple', 'Bm', 'Bp', 'Prot', 'cv']

        single_P1.iloc[:, 2:] = single_P1.iloc[:, 2:].mask(single_P1.iloc[:, 2:] == 0)
        single_P2.iloc[:, 2:] = single_P2.iloc[:, 2:].mask(single_P2.iloc[:, 2:] == 0)

        single_P1 = single_P1.dropna(axis=0)
        single_P2 = single_P2.dropna(axis=0)

        single_df = pd.concat([single_P1, single_P2])
        singletons = single_df.copy()

    else:
        singletons = pd.DataFrame(columns=['Round', 'Couple', 'Bm', 'Bp', 'Prot', 'cv'])
        duplicates = div_data.copy()

    return duplicates, singletons


def assign_paralogs(sim_data):
    """Function that generates lists of paralogs P1 (most expressed) and P2 (least expressed) from the
    simulation data

    sim_data = 'data' dataframe saved at the end of the simulation"""

    end_step = sim_data['Round'].unique().max()
    end_data = sim_data[sim_data['Round'] == end_step].copy()

    # Couples where the most expressed paralog (P1) is the one already identified as P1
    P1_left = np.where(end_data['Prot1'] >= end_data['Prot2'])

    # Couples where the most expressed paralog (P2) is the one that was identified as P2
    P1_right = np.where(end_data['Prot2'] >= end_data['Prot1'])

    return P1_left, P1_right, end_data


def divergence_panel(sim_data, n_couples, model_name, model_desc, loss=True):
    """Function to construct a 4-panel figure showing the divergence of all paralog couples in transcription rate,
    translation rate, protein abundance and expression noise (coefficient of variation) throughout the simulation.

    sim_data = 'data' dataframe saved at the end of the simulation
    n_couples = Number of paralog couples simulated
    model_name = Name of the model under which the current simulation has been done (to label the plots)
    model_desc = Longer description of the model used (to label the plots)
    loss = Boolean indicating whether gene loss was allowed during the simulation run used"""

    # All np.inf data is replaced by zeros
    sim_log2 = sim_data.copy()
    sim_log2.iloc[:, 2:] = sim_log2.iloc[:, 2:].where(sim_log2.iloc[:, 2:] != np.inf, other=0)

    # Then, all zeros are transformed into NaNs
    sim_log2.iloc[:, 2:] = sim_log2.iloc[:, 2:].where(sim_log2.iloc[:, 2:] != 0, other=np.NaN)

    # Log2 fold-changes are first calculated through time
    init_rates = sim_data[sim_data['Round'] == 0]
    n_rounds = sim_data['Round'].max() + 1
    init_tiled = np.tile(init_rates, (n_rounds, 1))

    sim_log2.iloc[:, 2:10] = sim_log2.iloc[:, 2:10] / init_tiled[:, 2:10]
    sim_log2.iloc[:, 2:10] = np.log2(sim_log2.iloc[:, 2:10])

    # Process data to separate singletons from remaining duplicates
    data_sep = process_data(sim_data, loss=loss)
    singletons = data_sep[1]

    # Most expressed paralogs (P1) are identified as being either on the left (identified as '1') or on the right of
    # the dataframe (identified as '2')
    paralogs = assign_paralogs(sim_data)
    P1_left = paralogs[0]
    P1_right = paralogs[1]

    # Data for P1 and P2 through time are aggregated
    means_P1 = pd.DataFrame(columns=['Round', 'Couple', 'Bm', 'Bp', 'Prot', 'cv'])
    means_P2 = means_P1.copy()

    for P1 in list(P1_left[0]):
        if P1 in singletons['Couple'].unique():
            singles = singletons[singletons['Couple'] == P1]
            couple_end = singles['Round'].min()
            subset = sim_log2[sim_log2['Couple'] == P1].copy()
            subset = subset[subset['Round'] < couple_end]

        else:
            subset = sim_log2[sim_log2['Couple'] == P1]

        P1_subset = subset.iloc[:, np.r_[0:4, 6, 8]].copy()
        P2_subset = subset.iloc[:, np.r_[0, 1, 4, 5, 7, 9]].copy()
        P1_subset.columns = ['Round', 'Couple', 'Bm', 'Bp', 'Prot', 'cv']
        P2_subset.columns = ['Round', 'Couple', 'Bm', 'Bp', 'Prot', 'cv']

        means_P1 = pd.concat([means_P1, P1_subset])
        means_P2 = pd.concat([means_P2, P2_subset])

    for P1 in list(P1_right[0]):
        if P1 in singletons['Couple'].unique():
            singles = singletons[singletons['Couple'] == P1]
            couple_end = singles['Round'].min()
            subset = sim_log2[sim_log2['Couple'] == P1].copy()
            subset = subset[subset['Round'] < couple_end]

        else:
            subset = sim_log2[sim_log2['Couple'] == P1]

        P1_subset = subset.iloc[:, np.r_[0, 1, 4, 5, 7, 9]].copy()
        P2_subset = subset.iloc[:, np.r_[0:4, 6, 8]].copy()
        P1_subset.columns = ['Round', 'Couple', 'Bm', 'Bp', 'Prot', 'cv']
        P2_subset.columns = ['Round', 'Couple', 'Bm', 'Bp', 'Prot', 'cv']

        means_P1 = pd.concat([means_P1, P1_subset])
        means_P2 = pd.concat([means_P2, P2_subset])

    # The means of both duplicates for each time point as well as the individual means of P1 and P2 are calculated
    means_both = pd.concat([means_P1, means_P2])
    means_both.iloc[:, 2:] = 2**means_both.iloc[:, 2:]
    means_both = means_both.groupby(by='Round', as_index=False).mean()
    means_both.iloc[:, 1:] = np.log2(means_both.iloc[0:, 1:])

    means_P1.iloc[:, 2:] = 2**means_P1.iloc[:, 2:]
    means_P1 = means_P1.groupby(by='Round', as_index=False).mean()
    means_P1.iloc[:, 1:] = np.log2(means_P1.iloc[0:, 1:])

    means_P2.iloc[:, 2:] = 2**means_P2.iloc[:, 2:]
    means_P2 = means_P2.groupby(by='Round', as_index=False).mean()
    means_P2.iloc[:, 1:] = np.log2(means_P2.iloc[0:, 1:])

    # The figure is constructed

    widths = 0.10

    if n_couples > 500:
        opacity = 0.5

    else:
        opacity = 0.9

    fig, axs = plt.subplots(2, 2, figsize=[16, 16])

    for P1 in list(P1_left[0]):
        # Here, paralog 1 (most expressed) is on the left of the df (identified as '1')
        subset = sim_log2[sim_log2['Couple'] == P1].copy()

        if P1 in singletons['Couple'].unique():
            singles = singletons[singletons['Couple'] == P1]
            couple_end = singles['Round'].min()
            subset = subset[subset['Round'] < couple_end]

        axs[0, 0].plot(subset['Round'], subset['Bm1'], c=cm.tab20.colors[3], alpha=opacity, linewidth=widths)
        axs[0, 0].plot(subset['Round'], subset['Bm2'], c=cm.tab20.colors[1], alpha=opacity, linewidth=widths)

        axs[0, 1].plot(subset['Round'], subset['Bp1'], c=cm.tab20.colors[3], alpha=opacity, linewidth=widths)
        axs[0, 1].plot(subset['Round'], subset['Bp2'], c=cm.tab20.colors[1], alpha=opacity, linewidth=widths)

        axs[1, 0].plot(subset['Round'], subset['Prot1'], c=cm.tab20.colors[3], alpha=opacity, linewidth=widths)
        axs[1, 0].plot(subset['Round'], subset['Prot2'], c=cm.tab20.colors[1], alpha=opacity, linewidth=widths)

        axs[1, 1].plot(subset['Round'], subset['cv1'], c=cm.tab20.colors[3], alpha=opacity, linewidth=widths)
        axs[1, 1].plot(subset['Round'], subset['cv2'], c=cm.tab20.colors[1], alpha=opacity, linewidth=widths)

    for P1 in list(P1_right[0]):
        # Here, paralog 1 (most expressed) is on the right of the df (identified as '2')
        subset = sim_log2[sim_log2['Couple'] == P1].copy()

        if P1 in singletons['Couple'].unique():
            singles = singletons[singletons['Couple'] == P1]
            couple_end = singles['Round'].min()
            subset = subset[subset['Round'] < couple_end]

        axs[0, 0].plot(subset['Round'], subset['Bm2'], c=cm.tab20.colors[3], alpha=opacity, linewidth=widths)
        axs[0, 0].plot(subset['Round'], subset['Bm1'], c=cm.tab20.colors[1], alpha=opacity, linewidth=widths)

        axs[0, 1].plot(subset['Round'], subset['Bp2'], c=cm.tab20.colors[3], alpha=opacity, linewidth=widths)
        axs[0, 1].plot(subset['Round'], subset['Bp1'], c=cm.tab20.colors[1], alpha=opacity, linewidth=widths)

        axs[1, 0].plot(subset['Round'], subset['Prot2'], c=cm.tab20.colors[3], alpha=opacity, linewidth=widths)
        axs[1, 0].plot(subset['Round'], subset['Prot1'], c=cm.tab20.colors[1], alpha=opacity, linewidth=widths)

        axs[1, 1].plot(subset['Round'], subset['cv2'], c=cm.tab20.colors[3], alpha=opacity, linewidth=widths)
        axs[1, 1].plot(subset['Round'], subset['cv1'], c=cm.tab20.colors[1], alpha=opacity, linewidth=widths)

    # Trajectories of singletons, if applicable
    if loss:

        # Dataframe to keep the log2 fold-changes of singleton genes
        single_folds = pd.DataFrame(columns=['Round', 'Couple', 'Bm', 'Bp', 'Prot', 'cv'])

        for gene in singletons['Couple'].unique():

            # The log2 fold-changes for singletons are obtained
            single_subset = singletons[singletons['Couple'] == gene]
            single_start = single_subset['Round'].min()

            subset = sim_log2[sim_log2['Couple'] == gene].copy()
            subset = subset[subset['Round'] >= single_start]

            # They are combined in a dataframe without the NaNs of the lost copies
            single_P1 = subset.iloc[:, np.r_[0:4, 6, 8]].copy()
            single_P2 = subset.iloc[:, np.r_[0:2, 4, 5, 7, 9]].copy()

            single_P1.columns = ['Round', 'Couple', 'Bm', 'Bp', 'Prot', 'cv']
            single_P2.columns = ['Round', 'Couple', 'Bm', 'Bp', 'Prot', 'cv']

            single_P1.iloc[:, 2:] = single_P1.iloc[:, 2:].mask(single_P1.iloc[:, 2:] == 0)
            single_P2.iloc[:, 2:] = single_P2.iloc[:, 2:].mask(single_P2.iloc[:, 2:] == 0)

            single_P1 = single_P1.dropna(axis=0)
            single_P2 = single_P2.dropna(axis=0)

            single_only = pd.concat([single_P1, single_P2])

            # Dataframe containing all singletons in log2 fold-changes
            single_folds = pd.concat([single_folds, single_only])

        # Plotting of the individual trajectories
        for gene in single_folds['Couple'].unique():

            folds_sub = single_folds[single_folds['Couple'] == gene]

            axs[0, 0].plot(folds_sub['Round'], folds_sub['Bm'], c=cm.tab20.colors[7], alpha=opacity,
                           linewidth=widths)
            axs[0, 1].plot(folds_sub['Round'], folds_sub['Bp'], c=cm.tab20.colors[7], alpha=opacity,
                           linewidth=widths)
            axs[1, 0].plot(folds_sub['Round'], folds_sub['Prot'], c=cm.tab20.colors[7], alpha=opacity,
                           linewidth=widths)
            axs[1, 1].plot(folds_sub['Round'], folds_sub['cv'], c=cm.tab20.colors[7], alpha=opacity,
                           linewidth=widths)

    # P1 means
    axs[0, 0].plot(means_P1['Round'], means_P1['Bm'], c=cm.tab10.colors[1], label='P1')
    axs[0, 1].plot(means_P1['Round'], means_P1['Bp'], c=cm.tab10.colors[1], label='P1')
    axs[1, 0].plot(means_P1['Round'], means_P1['Prot'], c=cm.tab10.colors[1], label='P1')
    axs[1, 1].plot(means_P1['Round'], means_P1['cv'], c=cm.tab10.colors[1], label='P1')

    # P2 means
    axs[0, 0].plot(means_P2['Round'], means_P2['Bm'], c=cm.tab10.colors[0], label='P2')
    axs[0, 1].plot(means_P2['Round'], means_P2['Bp'], c=cm.tab10.colors[0], label='P2')
    axs[1, 0].plot(means_P2['Round'], means_P2['Prot'], c=cm.tab10.colors[0], label='P2')
    axs[1, 1].plot(means_P2['Round'], means_P2['cv'], c=cm.tab10.colors[0], label='P2')

    # Global means
    axs[0, 0].plot(means_both['Round'], means_both['Bm'], c='k', label='Duplicates', alpha=0.65)
    axs[0, 1].plot(means_both['Round'], means_both['Bp'], c='k', label='Duplicates', alpha=0.65)
    axs[1, 0].plot(means_both['Round'], means_both['Prot'], c='k', label='Duplicates', alpha=0.65)
    axs[1, 1].plot(means_both['Round'], means_both['cv'], c='k', label='Duplicates', alpha=0.65)

    # Mean of singletons, if applicable
    if loss:
        means_single = single_folds.copy()
        means_single.iloc[:, 1:] = 2 ** means_single.iloc[:, 1:]
        means_single = means_single.groupby(by='Round', as_index=False).mean()
        means_single.iloc[:, 1:] = np.log2(means_single.iloc[0:, 1:])

        axs[0, 0].plot(means_single['Round'], means_single['Bm'], c=cm.tab10.colors[3], label='Singletons')
        axs[0, 1].plot(means_single['Round'], means_single['Bp'], c=cm.tab10.colors[3], label='Singletons')
        axs[1, 0].plot(means_single['Round'], means_single['Prot'], c=cm.tab10.colors[3], label='Singletons')
        axs[1, 1].plot(means_single['Round'], means_single['cv'], c=cm.tab10.colors[3], label='Singletons')

    # Labels
    for ax in [axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]:
        ax.legend()
        ax.set_xlabel('Mutation Round')
        ax.set_ylabel('Log2 fold-change')

    axs[0, 0].set_title('Transcription rate')
    axs[0, 1].set_title('Translation rate')
    axs[1, 0].set_title('Protein abundance')
    axs[1, 1].set_title('Expression noise (cv)')
    fig.suptitle(f'Individual trajectories under the {model_name} model\n({model_desc})', fontsize=14)

    return fig


def trajectories_space(sim_data, n_couples, total_bm, min_bm, max_bm, min_bp, max_bp, loss=False,
                       model='mixed', boundary=2.0544):
    """Function to plot the trajectories (initial position -> final position) of each replicate duplicate couple
    in the expression space (Log10 transcription and translation rates)

    sim_data = 'data' dataframe saved at the end of the simulation
    n_couples = Number of replicate couples in the current simulation run
    total_bm = Ratio of total transcriptional flux after duplication over the initial transcriptional flux. If it is
               not 2, it means that the initial position of the duplicates is not the same as the position of the
               ancestral singleton.
    loss = Boolean indicating whether gene loss was allowed during the simulation run used
    boundary = Bp/Bm ratio defining the depleted region of the expression space, as reported by Hausser et al., 2019"""

    # Initial transcription and translation rates
    post_bm = list(sim_data.iloc[0:n_couples, 2])
    post_bp = list(sim_data.iloc[0:n_couples, 3])

    # Process data to separate singletons from remaining duplicates
    data_sep = process_data(sim_data, loss=loss)
    dupli_df = data_sep[0]
    singletons = data_sep[1]

    # Most expressed paralogs (P1) are identified as being either on the left (identified as '1') or on the right of
    # the dataframe (identified as '2')
    paralogs = assign_paralogs(sim_data)
    P1_left = paralogs[0]
    P1_right = paralogs[1]

    end_step = sim_data['Round'].unique().max()
    end_dupli = dupli_df[dupli_df['Round'] == end_step].copy()

    # The figure is constructed

    widths = 0.75

    if n_couples > 500:
        opacity = 0.75

    else:
        opacity = 1.0

    fig = plt.figure(figsize=(9, 9))

    # All trajectories for paralogs P1 (most expressed) and P2 (least expressed) of duplicate couples
    if loss:
        end_single = singletons[singletons['Round'] == end_step].copy()
        plt.plot(0, 0, c=cm.tab20.colors[7], label='Singletons')

    else:
        end_single = None

    for couple in range(n_couples):
        if couple in end_dupli['Couple'].unique():

            if couple in P1_left[0]:
                bm_traj_P1 = np.log10(np.array([post_bm[couple], end_dupli[end_dupli['Couple'] == couple]['Bm1'].values[0]]))
                bp_traj_P1 = np.log10(np.array([post_bp[couple], end_dupli[end_dupli['Couple'] == couple]['Bp1'].values[0]]))
                bm_traj_P2 = np.log10(np.array([post_bm[couple], end_dupli[end_dupli['Couple'] == couple]['Bm2'].values[0]]))
                bp_traj_P2 = np.log10(np.array([post_bp[couple], end_dupli[end_dupli['Couple'] == couple]['Bp2'].values[0]]))

                plt.plot(bm_traj_P1, bp_traj_P1, c=cm.tab20.colors[3], alpha=opacity, linewidth=widths)
                plt.plot(bm_traj_P2, bp_traj_P2, c=cm.tab20.colors[1], alpha=opacity, linewidth=widths)

            elif couple in P1_right[0]:
                bm_traj_P1 = np.log10(np.array([post_bm[couple], end_dupli[end_dupli['Couple'] == couple]['Bm2'].values[0]]))
                bp_traj_P1 = np.log10(np.array([post_bp[couple], end_dupli[end_dupli['Couple'] == couple]['Bp2'].values[0]]))
                bm_traj_P2 = np.log10(np.array([post_bm[couple], end_dupli[end_dupli['Couple'] == couple]['Bm1'].values[0]]))
                bp_traj_P2 = np.log10(np.array([post_bp[couple], end_dupli[end_dupli['Couple'] == couple]['Bp1'].values[0]]))

                plt.plot(bm_traj_P1, bp_traj_P1, c=cm.tab20.colors[3], alpha=opacity, linewidth=widths)
                plt.plot(bm_traj_P2, bp_traj_P2, c=cm.tab20.colors[1], alpha=opacity, linewidth=widths)

        elif couple in end_single['Couple'].unique():
            bm_traj = np.log10(np.array([post_bm[couple], end_single[end_single['Couple'] == couple]['Bm'].values[0]]))
            bp_traj = np.log10(np.array([post_bp[couple], end_single[end_single['Couple'] == couple]['Bp'].values[0]]))

            plt.plot(bm_traj, bp_traj, c=cm.tab20.colors[7], alpha=opacity, linewidth=widths)
            plt.plot(bm_traj, bp_traj, c=cm.tab20.colors[7], alpha=opacity, linewidth=widths)

    plt.plot(0, 0, c=cm.tab20.colors[3], label='P1')
    plt.plot(0, 0, c=cm.tab20.colors[1], label='P2')
    plt.legend()

    # Points representing the initial locations of the replicates
    for couple in range(n_couples):
        plt.scatter(math.log10(post_bm[couple]), math.log10(post_bp[couple]), color='black', s=1,
                    label='Post-duplication')

    # Addition of the boundary of the depleted region of the expression space, as reported by Hausser et al., 2019
    log10_k = math.log10(boundary)
    bm_interval = np.linspace(-1, 4, 50)
    plt.plot(bm_interval, bm_interval + log10_k, '-k')

    # Addition of the boundaries of minimal and maximal transcription and translation rates
    if min_bm != 0:
        plt.axvline(x=np.log10(min_bm), c='c', linestyle='--')
    plt.axvline(x=np.log10(max_bm), c='c', linestyle='--')

    if min_bp != 0:
        plt.axhline(y=np.log10(min_bp), c='k', linestyle='--')
    plt.axhline(y=np.log10(max_bp), c='k', linestyle='--')

    # The percentage of genes in the depleted region is assessed
    P1_out = np.where(end_dupli['Bp1']/end_dupli['Bm1'] < boundary, 1, 0)
    P2_out = np.where(end_dupli['Bp2']/end_dupli['Bm2'] < boundary, 1, 0)

    if loss:
        single_out = np.where(end_single['Bp']/end_single['Bm'] < boundary, 1, 0)
        all_out = np.sum(P1_out) + np.sum(P2_out) + np.sum(single_out)

    else:
        all_out = np.sum(P1_out) + np.sum(P2_out)

    percent_out = (all_out/(2 * n_couples)) * 100

    plt.annotate(f'{round(percent_out, 2)} % of\nall genes', (3, 1.5), fontsize=12, bbox=dict(boxstyle='round', fc='w'))

    plt.xlim(left=-1, right=4)
    plt.xlabel('Log10 transcription rate [mRNA/h]')

    plt.ylim(bottom=0)
    plt.ylabel('Log10 translation rate [protein / (mRNA*h)')

    if total_bm != 2.0:
        plt.title(f'Evolutionary trajectories in the expression space for the {model} model\n(The post-duplication '
                  f'state does not reflect the ancestral singleton)')

    else:
        plt.title(f'Evolutionary trajectories in the expression space for the {model} model')

    return fig
