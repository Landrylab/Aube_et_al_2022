# This script combines multiple iterations of the simulation to produce a dataframe of mean KS statistics.
import os
import glob
import pandas as pd

iters = ['effect_2', 'effect_16', 'effect_186', 'effect_225']

effects = {'effect_2': 2.0, 'effect_16': 1.6, 'effect_186': 1.86, 'effect_225': 2.25}

KS_all = pd.DataFrame(columns=['Model', 'Comparison', 'Property', 'KS stats', 'KS p-value', 'Run', 'Dupli_effect'])

for directory in iters:
    os.chdir(f'F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2/Dupli_effect/{directory}')

    df_KS = pd.concat(map(pd.read_csv, glob.glob('KS_tests*.csv')))
    df_KS['Dupli_effect'] = effects[directory]

    KS_all = pd.concat([KS_all, df_KS])

os.chdir('F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2/Dupli_effect')

# 'iter{number}' are deleted from the Run tags
KS_all = KS_all.reset_index(drop=True)
for row in range(KS_all.shape[0]):
    KS_all.at[row, 'Run'] = KS_all.at[row, 'Run'][:-6]

# Mean KS statistics are computed
KS_mean = KS_all.groupby(by=['Model', 'Comparison', 'Run', 'Dupli_effect'], as_index=False).mean()
KS_mean = KS_mean[KS_mean['Comparison'] == 'All duplicates'].reset_index(drop=True)
KS_mean = KS_mean.drop('Comparison', axis=1)

KS_mean = KS_mean.astype({'Dupli_effect':'float64'})

KS_mean.to_csv('KS_dupli_effect.csv', index=False)
