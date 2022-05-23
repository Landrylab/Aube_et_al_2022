import os
import glob
import numpy as np
import pandas as pd
import scipy.stats as stats

# 1) For Mood's p-values
iters = ['effect_2', 'effect_16', 'effect_186', 'effect_225']

dupli_effect = {'effect_2': 2, 'effect_16': 1.6, 'effect_186': 1.86, 'effect_225': 2.25}

folds_model = pd.DataFrame(columns=['Model', 'Parameter', 'Log2 fold-change', 'Dupli_effect', 'Run'])
folds_all = folds_model.copy()

for directory in iters:
    os.chdir(f'F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2'
             f'/Dupli_effect/{directory}')

    folds_direct = folds_model.copy()

    for bias in ['Bm05', 'Bm1', 'Bm2', 'Bm3', 'Bm4', 'Bm5', 'Bm6', 'Bm7', 'Bm8', 'Bm9', 'Bm10']:

        for run in [1, 2, 3]:

            if os.path.exists(f'Fold_changes_all_{bias}_iter{run}.csv'):
                df = pd.read_csv(f'Fold_changes_all_{bias}_iter{run}.csv')
                df['Run'] = f'{bias}_iter{run}'

                folds_direct = pd.concat([folds_direct, df])

    folds_direct['Dupli_effect'] = dupli_effect[directory]

    folds_all = pd.concat([folds_all, folds_direct])

os.chdir('F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2/Dupli_effect')

# Simulation data and yeast data are separated:
yeast_dist = folds_all[folds_all['Model'] == 'WGD couples']
# Only one of the identical series of values (the same values had been appended to each df) is kept
yeast_dist = yeast_dist[(yeast_dist['Run'] == yeast_dist['Run'].unique()[0]) &
                        (yeast_dist['Dupli_effect'] == yeast_dist['Dupli_effect'].unique()[0])]

sim_dist = folds_all[(folds_all['Model'] == 'Mixed') | (folds_all['Model'] == 'ADS') |
                     (folds_all['Model'] == 'Minimal') | (folds_all['Model'] == 'No Cost')]

# Mood's median test is performed and p-values as well as test statistics are kept
Moods_res = pd.DataFrame(columns=['Model', 'Dupli_effect', 'Run', 'Parameter', 'p-value', 'Statistics'])

yeast_bm = yeast_dist[yeast_dist['Parameter'] == 'Transcription rate']
yeast_bp = yeast_dist[yeast_dist['Parameter'] == 'Translation rate']
yeast_prot = yeast_dist[yeast_dist['Parameter'] == 'Protein abundance']
yeast_ratio = yeast_dist[yeast_dist['Parameter'] == 'Divergence ratio']

row = 0

for corr in sim_dist['Dupli_effect'].unique():
    corr_sub = sim_dist[sim_dist['Dupli_effect'] == corr]

    for run in corr_sub['Run'].unique():
        run_sub = corr_sub[corr_sub['Run'] == run]

        for model in ['Mixed', 'ADS', 'Minimal', 'No Cost']:
            model_sub = run_sub[run_sub['Model'] == model]

            bm = model_sub[model_sub['Parameter'] == 'Transcription rate']
            bp = model_sub[model_sub['Parameter'] == 'Translation rate']
            prot = model_sub[model_sub['Parameter'] == 'Protein abundance']
            ratio = model_sub[model_sub['Parameter'] == 'Divergence ratio']

            mood_bm = stats.median_test(bm['Log2 fold-change'], yeast_bm['Log2 fold-change'])
            Moods_res.loc[row] = [f'{model}', f'{corr}', f'{run}', 'Transcription rate', mood_bm[1], mood_bm[0]]

            row += 1

            mood_bp = stats.median_test(bp['Log2 fold-change'], yeast_bp['Log2 fold-change'])
            Moods_res.loc[row] = [f'{model}', f'{corr}', f'{run}', 'Translation rate', mood_bp[1], mood_bp[0]]

            row += 1

            mood_prot = stats.median_test(prot['Log2 fold-change'], yeast_prot['Log2 fold-change'])
            Moods_res.loc[row] = [f'{model}', f'{corr}', f'{run}', 'Protein abundance', mood_prot[1], mood_prot[0]]

            row += 1

            mood_ratio = stats.median_test(ratio['Log2 fold-change'], yeast_ratio['Log2 fold-change'])
            Moods_res.loc[row] = [f'{model}', f'{corr}', f'{run}', 'Divergence ratio', mood_ratio[1], mood_ratio[0]]

            row += 1

# A new iter column is added, and the 'iter' mention
# is deleted from the 'Run column
Moods_res['Iter'] = 'problem'

for row in range(Moods_res.shape[0]):
    Moods_res.at[row, 'Iter'] = Moods_res.at[row, 'Run'][-5:]
    Moods_res.at[row, 'Run'] = Moods_res.at[row, 'Run'][:-6]

# The dataframe of Mood's p-values is exported
Moods_res.to_csv('Mood_results.csv', index=False)

# 2) For KS statistics
KS_all = pd.DataFrame(columns=['Model', 'Comparison', 'Property', 'KS stats', 'KS p-value', 'Run', 'Dupli_effect'])

for directory in iters:
    os.chdir(f'F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2'
             f'/Dupli_effect/{directory}')

    df_KS = pd.concat(map(pd.read_csv, glob.glob('KS_tests*.csv')))
    df_KS['Dupli_effect'] = dupli_effect[directory]

    KS_all = pd.concat([KS_all, df_KS])

os.chdir('F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2/Dupli_effect')

# 'iter{number}' are deleted from the Run tags
KS_all = KS_all.reset_index(drop=True)
for row in range(KS_all.shape[0]):
    KS_all.at[row, 'Run'] = KS_all.at[row, 'Run'][:-6]

# Mean KS statistics are computed
KS_bias = KS_all.groupby(by=['Model', 'Comparison', 'Property', 'Run', 'Dupli_effect'], as_index=False).mean()

# And KS statistics are multiplied by -1
KS_bias['KS stats'] = KS_bias['KS stats'] * -1

KS_bias = KS_bias.astype({'Dupli_effect':'float64'})

# The dataframe of Ks statistics is saved
KS_bias.to_csv('KS_dupli.csv', index=False)

# 3) For correlation coefficients

# Data are first imported iteratively from all simulation runs
data_all = pd.DataFrame(columns=['Round', 'Couple', 'Bm1', 'Bp1', 'Bm2', 'Bp2', 'Prot1', 'Prot2', 'Run', 'Model',
                                 'Dupli_effect'])
data_model = pd.DataFrame(columns=['Round', 'Couple', 'Bm1', 'Bp1', 'Bm2', 'Bp2', 'Prot1', 'Prot2', 'Run'])

for folder in iters:
    os.chdir(f'F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2/'
             f'Dupli_effect/{folder}')

    data_ADS = data_model.copy()
    files_ADS = glob.glob('data_ADS*.csv')

    for name in files_ADS:
        df = pd.read_csv(name)
        df['Run'] = name[9:-4]

        df = df[df['Round'] != 0]
        df = df.drop(labels=['cv1', 'cv2', 'Exp_cost'], axis=1)

        data_ADS = pd.concat([data_ADS, df])

    data_ADS['Model'] = 'ADS'
    data_ADS['Dupli_effect'] = dupli_effect[folder]
    data_ADS = data_ADS.reset_index(drop=True)

    # Same for the mixed data
    data_mixed = data_model.copy()
    files_mixed = glob.glob('data_Mixed*.csv')

    for name in files_mixed:
        df = pd.read_csv(name)
        df['Run'] = name[11:-4]

        df = df[df['Round'] != 0]
        df = df.drop(labels=['cv1', 'cv2', 'Exp_cost'], axis=1)

        data_mixed = pd.concat([data_mixed, df])

    data_mixed['Model'] = 'Mixed'
    data_mixed['Dupli_effect'] = dupli_effect[folder]
    data_mixed = data_mixed.reset_index(drop=True)

    # Same for minimal data
    data_minimal = data_model.copy()
    files_minimal = glob.glob('data_minimal*.csv')

    for name in files_minimal:
        df = pd.read_csv(name)
        df['Run'] = name[13:-4]

        df = df[df['Round'] != 0]
        df = df.drop(labels=['cv1', 'cv2', 'Exp_cost'], axis=1)

        data_minimal = pd.concat([data_minimal, df])

    data_minimal['Model'] = 'Minimal'
    data_minimal['Dupli_effect'] = dupli_effect[folder]
    data_minimal = data_minimal.reset_index(drop=True)

    # All data is combined
    data_all = pd.concat([data_all, data_ADS, data_mixed, data_minimal])

os.chdir('F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2/Dupli_effect')

# All transcriptional and translational log2 fold-changes are first calculated in the same
# direction (most expressed / least expressed), according to final protein abundance
data_all = data_all.reset_index(drop=True)

P1_higher = np.where(data_all['Prot1'] >= data_all['Prot2'], 1, 0)
P2_higher = np.where(data_all['Prot1'] < data_all['Prot2'], 1, 0)

Bm1_higher = np.where(data_all['Bm1'] >= data_all['Bm2'], 1, 0)
Bm2_higher = np.where(data_all['Bm1'] < data_all['Bm2'], 1, 0)
Bp1_higher = np.where(data_all['Bp1'] >= data_all['Bp2'], 1, 0)
Bp2_higher = np.where(data_all['Bp1'] < data_all['Bp2'], 1, 0)

data_P1 = data_all.copy()
data_P2 = data_all.copy()

data_P1['Logfold_Bm'] = data_P1['Bm1'] / data_P1['Bm2']
data_P1['Logfold_Bp'] = data_P1['Bp1'] / data_P1['Bp2']
data_P1['Logfold_pEst'] = data_P1['Prot1'] / data_P1['Prot2']

data_P2['Logfold_Bm'] = data_P2['Bm2'] / data_P2['Bm1']
data_P2['Logfold_Bp'] = data_P2['Bp2'] / data_P2['Bp1']
data_P2['Logfold_pEst'] = data_P2['Prot2'] / data_P2['Prot1']

data_P1.iloc[:, 12] = data_P1.iloc[:, 12] * Bm1_higher
data_P1.iloc[:, 13] = data_P1.iloc[:, 13] * Bp1_higher
data_P1.iloc[:, 14] = data_P1.iloc[:, 14] * P1_higher

data_P2.iloc[:, 12] = data_P2.iloc[:, 12] * Bm2_higher
data_P2.iloc[:, 13] = data_P2.iloc[:, 13] * Bp2_higher
data_P2.iloc[:, 14] = data_P2.iloc[:, 14] * P2_higher

# All data (P1/P2 and P2/P1) put together, then converted to log2 fold-changes
log_both = data_P1.copy()
log_both.iloc[:, 12:15] = (log_both.iloc[:, 12:15] + data_P2.iloc[:, 12:15])
log_both['Logfold_Bm'] = np.log2(log_both['Logfold_Bm'])
log_both['Logfold_Bp'] = np.log2(log_both['Logfold_Bp'])
log_both['Logfold_pEst'] = np.log2(log_both['Logfold_pEst'])
log_both['Div_ratio'] = log_both['Logfold_Bm'] - log_both['Logfold_Bp']

# The desired correlations can be computed
corr_model = pd.DataFrame(columns=['Model', 'Run', 'Iter', 'Dupli_effect', 'pearson', 'r_pval', 'spearman', 'rho_pval'])
corr_divergence = corr_model.copy()
corr_sign = corr_model.copy()
corr_ratio = corr_model.copy()

for model in ['ADS', 'Mixed', 'Minimal', 'No Cost']:
    model_sub = log_both[log_both['Model'] == model]

    for corr in model_sub['Dupli_effect'].unique():
        corr_sub = model_sub[model_sub['Dupli_effect'] == corr]

        for run in corr_sub['Run'].unique():
            run_sub = corr_sub[corr_sub['Run'] == run]

            # First, divergence correlations
            corr_data_div = corr_model.copy()
            corr_data_div.at[0, 'Model'] = model
            corr_data_div.at[0, 'Iter'] = run[-5:]
            corr_data_div.at[0, 'Run'] = run[:-6]
            corr_data_div.at[0, 'Dupli_effect'] = corr

            pearson_div = stats.pearsonr(run_sub['Logfold_Bm'], run_sub['Logfold_Bp'])
            spearman_div = stats.spearmanr(run_sub['Logfold_Bm'], run_sub['Logfold_Bp'])

            corr_data_div.at[0, 'pearson'] = pearson_div[0]
            corr_data_div.at[0, 'r_pval'] = pearson_div[1]
            corr_data_div.at[0, 'spearman'] = spearman_div[0]
            corr_data_div.at[0, 'rho_pval'] = spearman_div[1]

            corr_divergence = pd.concat([corr_divergence, corr_data_div]).reset_index(drop=True)

            # Second, for correlations of signed divergence
            corr_data_sign = corr_model.copy()
            corr_data_sign.at[0, 'Model'] = model
            corr_data_sign.at[0, 'Iter'] = run[-5:]
            corr_data_sign.at[0, 'Run'] = run[:-6]
            corr_data_sign.at[0, 'Dupli_effect'] = corr

            Bm_signed = np.log2(run_sub['Bm1'] / run_sub['Bm2'])
            Bp_signed = np.log2(run_sub['Bp1'] / run_sub['Bp2'])

            pearson_sign = stats.pearsonr(Bm_signed, Bp_signed)
            spearman_sign = stats.spearmanr(Bm_signed, Bp_signed)

            corr_data_sign.at[0, 'pearson'] = pearson_sign[0]
            corr_data_sign.at[0, 'r_pval'] = pearson_sign[1]
            corr_data_sign.at[0, 'spearman'] = spearman_sign[0]
            corr_data_sign.at[0, 'rho_pval'] = spearman_sign[1]

            corr_sign = pd.concat([corr_sign, corr_data_sign]).reset_index(drop=True)

            # Third, for correlations with the divergence ratio
            corr_data_ratio = corr_model.copy()
            corr_data_ratio.at[0, 'Model'] = model
            corr_data_ratio.at[0, 'Iter'] = run[-5:]
            corr_data_ratio.at[0, 'Run'] = run[:-6]
            corr_data_ratio.at[0, 'Dupli_effect'] = corr

            pearson_ratio = stats.pearsonr(run_sub['Logfold_pEst'], run_sub['Div_ratio'])
            spearman_ratio = stats.spearmanr(run_sub['Logfold_pEst'], run_sub['Div_ratio'])

            corr_data_ratio.at[0, 'pearson'] = pearson_ratio[0]
            corr_data_ratio.at[0, 'r_pval'] = pearson_ratio[1]
            corr_data_ratio.at[0, 'spearman'] = spearman_ratio[0]
            corr_data_ratio.at[0, 'rho_pval'] = spearman_ratio[1]

            corr_ratio = pd.concat([corr_ratio, corr_data_ratio]).reset_index(drop=True)

# The dataframe of correlations is exported
corr_divergence.to_csv('Corr_dupli.csv', index=False)
corr_sign.to_csv('Corr_sign_dupli.csv', index=False)
corr_ratio.to_csv('Corr_ratio_dupli.csv', index=False)
