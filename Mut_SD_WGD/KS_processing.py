# Short script to process the KS statistics for the later construction of Supplementary Fig. 3
import os
import glob
import pandas as pd

# Mapping of directory names and sigma_mut
iters = ['1_percent', '2_5_percent', '5_percent', '7_5_percent', '10_percent', '12_5_percent', '15_percent',
         '20_percent', '25_percent', '35_percent']

sd_vals = {'1_percent': 0.01, '2_5_percent': 0.025, '5_percent': 0.05, '7_5_percent': 0.075, '10_percent': 0.1,
           '12_5_percent': 0.125, '15_percent': 0.15, '20_percent': 0.2, '25_percent': 0.25, '35_percent': 0.35}

# Combining all data in one dataframe
KS_all = pd.DataFrame(columns=['Model', 'Comparison', 'Property', 'KS stats', 'KS p-value', 'Run', 'Mut_sd'])

for directory in iters:
    os.chdir(f'F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/'
             f'Final_simulations_v2/Mut_SD_WGD/{directory}')

    df_KS = pd.concat(map(pd.read_csv, glob.glob('KS_tests*.csv')))
    df_KS['Mut_sd'] = sd_vals[directory]

    KS_all = pd.concat([KS_all, df_KS])

os.chdir(f'F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2/Mut_SD_WGD')

# 'iter{number}' are deleted from the Run tags
KS_all = KS_all.reset_index(drop=True)
for row in range(KS_all.shape[0]):
    KS_all.at[row, 'Run'] = KS_all.at[row, 'Run'][:-6]

# Only the comparison with WGD duplicate couples is kept
KS_all = KS_all[KS_all['Comparison'] == 'WGD'].copy().reset_index(drop=True)

# Generation of a KS_modified.csv dataframe to be used for the Supp figure 4
KS_sub = KS_all[KS_all['Mut_sd'] == 0.100].copy().reset_index(drop=True)
mut_ratio = {'Bm05': '1/2', 'Bm1': '1', 'Bm2': '2', 'Bm3': '3', 'Bm4': '4', 'Bm5': '5', 'Bm6': '6', 'Bm7': '7',
             'Bm8': '8', 'Bm9': '9', 'Bm10': '10'}

for row in range(KS_sub.shape[0]):
    KS_sub.at[row, 'Run'] = mut_ratio[KS_sub.at[row, 'Run']]

KS_sub = KS_sub.rename(columns={'Run': 'Ratio'})
KS_sub.to_csv('KS_modified_WGD.csv')

# Mean KS statistics are computed
KS_all = KS_all.drop('KS p-value', axis=1)
KS_means = KS_all.groupby(by=['Model', 'Comparison', 'Run', 'Mut_sd'], as_index=False).mean()

# The dataframe is exported as csv to be used in the Colaboratory notebook where the Supp. figure is generated
KS_means.to_csv('KS_means_WGD.csv')
