# This script combines multiple iterations of the simulation to produce heatmaps of Mood's p-values and KS
# statistics according to two parameters.
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
import scipy.stats as stats
from matplotlib.backends.backend_pdf import PdfPages

# 1) For Mood's p-values
iters = ['corr_0', 'corr_01', 'corr_01_neg', 'corr_02', 'corr_02_neg', 'corr_03', 'corr_03_neg',
         'corr_04', 'corr_04_neg', 'corr_05', 'corr_06']

corr_vals = {'corr_0': 0, 'corr_01': 0.1, 'corr_01_neg': -0.1, 'corr_02': 0.2, 'corr_02_neg': -0.2,
             'corr_03': 0.3, 'corr_03_neg': -0.3, 'corr_04': 0.4, 'corr_04_neg': -0.4, 'corr_05': 0.5,
             'corr_06': 0.6, 'corr_07': 0.7}

folds_model = pd.DataFrame(columns=['Model', 'Parameter', 'Log2 fold-change', 'Mut_corr', 'Run'])
folds_all = folds_model.copy()


for directory in iters:
    os.chdir(f'F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2/'
             f'Mut_correlations_WGD_1e5/{directory}')

    folds_direct = folds_model.copy()

    for bias in ['Bm05', 'Bm1', 'Bm2', 'Bm3', 'Bm4', 'Bm5', 'Bm6', 'Bm7', 'Bm8', 'Bm9', 'Bm10']:

        for run in [1, 2, 3]:

            if os.path.exists(f'Fold_changes_all_{bias}_iter{run}.csv'):
                df = pd.read_csv(f'Fold_changes_all_{bias}_iter{run}.csv')
                df['Run'] = f'{bias}_iter{run}'

                folds_direct = pd.concat([folds_direct, df])

    folds_direct['Mut_corr'] = corr_vals[directory]

    folds_all = pd.concat([folds_all, folds_direct])

os.chdir('F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2/'
         'Mut_correlations_WGD_1e5')

# Simulation data and yeast data are separated:
yeast_dist = folds_all[folds_all['Model'] == 'Yeast Data']
# Only one of the identical series of values (the same values had been appended to each df) is kept
yeast_dist = yeast_dist[(yeast_dist['Run'] == yeast_dist['Run'].unique()[0]) &
                        (yeast_dist['Mut_corr'] == yeast_dist['Mut_corr'].unique()[0])]

sim_dist = folds_all[(folds_all['Model'] == 'Mixed') | (folds_all['Model'] == 'ADS') | (folds_all['Model'] == 'Minimal')]

# Mood's median test is performed and p-values as well as test statistics are kept
Moods_res = pd.DataFrame(columns=['Model', 'Mut_corr', 'Run', 'Parameter', 'p-value', 'Statistics'])

yeast_bm = yeast_dist[yeast_dist['Parameter'] == 'Transcription rate']
yeast_bp = yeast_dist[yeast_dist['Parameter'] == 'Translation rate']
yeast_prot = yeast_dist[yeast_dist['Parameter'] == 'Protein abundance']
yeast_ratio = yeast_dist[yeast_dist['Parameter'] == 'Divergence ratio']

row = 0

for corr in sim_dist['Mut_corr'].unique():
    corr_sub = sim_dist[sim_dist['Mut_corr'] == corr]

    for run in corr_sub['Run'].unique():
        run_sub = corr_sub[corr_sub['Run'] == run]

        for model in ['Mixed', 'ADS', 'Minimal']:
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

# The 'iter' mention is deleted from the 'Run column
for row in range(Moods_res.shape[0]):

    Moods_res.at[row, 'Run'] = Moods_res.at[row, 'Run'][:-6]

# The matrices are prepared for the construction of the heatmaps
Moods_means = Moods_res.groupby(by=['Model', 'Mut_corr', 'Run', 'Parameter'], as_index=False).mean()
Moods_means['p-value'] = np.log10(Moods_means['p-value'])
Moods_means = Moods_means.astype({'Mut_corr': 'float64'})

bm_moods = Moods_means[Moods_means['Parameter'] == 'Transcription rate']
bp_moods = Moods_means[Moods_means['Parameter'] == 'Translation rate']
prot_moods = Moods_means[Moods_means['Parameter'] == 'Protein abundance']
ratio_moods = Moods_means[Moods_means['Parameter'] == 'Divergence ratio']

bm_mixed = bm_moods[bm_moods['Model'] == 'Mixed'].pivot(index='Mut_corr', columns='Run', values='p-value')
bp_mixed = bp_moods[bp_moods['Model'] == 'Mixed'].pivot(index='Mut_corr', columns='Run', values='p-value')
prot_mixed = prot_moods[prot_moods['Model'] == 'Mixed'].pivot(index='Mut_corr', columns='Run', values='p-value')
ratio_mixed = ratio_moods[ratio_moods['Model'] == 'Mixed'].pivot(index='Mut_corr', columns='Run', values='p-value')

bm_ADS = bm_moods[bm_moods['Model'] == 'ADS'].pivot(index='Mut_corr', columns='Run', values='p-value')
bp_ADS = bp_moods[bp_moods['Model'] == 'ADS'].pivot(index='Mut_corr', columns='Run', values='p-value')
prot_ADS = prot_moods[prot_moods['Model'] == 'ADS'].pivot(index='Mut_corr', columns='Run', values='p-value')
ratio_ADS = ratio_moods[ratio_moods['Model'] == 'ADS'].pivot(index='Mut_corr', columns='Run', values='p-value')

bm_min = bm_moods[bm_moods['Model'] == 'Minimal'].pivot(index='Mut_corr', columns='Run', values='p-value')
bp_min = bp_moods[bp_moods['Model'] == 'Minimal'].pivot(index='Mut_corr', columns='Run', values='p-value')
prot_min = prot_moods[prot_moods['Model'] == 'Minimal'].pivot(index='Mut_corr', columns='Run', values='p-value')
ratio_min = ratio_moods[ratio_moods['Model'] == 'Minimal'].pivot(index='Mut_corr', columns='Run', values='p-value')

# Column names are modified
dfs = [bm_mixed, bp_mixed, prot_mixed, ratio_mixed, bm_ADS, bp_ADS, prot_ADS, ratio_ADS, bm_min, bp_min, prot_min,
       ratio_min]
for df in dfs:
    df.columns = ['1/2', '1', '10', '2', '3', '4', '5', '6', '7', '8', '9']

# For reordering the matrix
columns = ['1/2', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

# Generation of the figure
fig, axs = plt.subplots(4, 3, figsize=[48, 48])

fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.93, 0.20, 0.02, 0.6])

sns.heatmap(bm_mixed.reindex(labels=columns, axis='columns'), vmin=-12, vmax=np.max(Moods_means['p-value']),
            cmap='bwr_r', linewidths=0.35, linecolor='grey', center=np.log10(0.05), ax=axs[0,0], cbar_ax=cbar_ax)
sns.heatmap(bm_ADS.reindex(labels=columns, axis='columns'), vmin=-12, vmax=np.max(Moods_means['p-value']),
            cmap='bwr_r', linewidths=0.35, linecolor='grey', center=np.log10(0.05), ax=axs[0,1], cbar=False)
sns.heatmap(bm_min.reindex(labels=columns, axis='columns'), vmin=-12, vmax=np.max(Moods_means['p-value']),
            cmap='bwr_r', linewidths=0.35, linecolor='grey', center=np.log10(0.05), ax=axs[0,2], cbar=False)

sns.heatmap(bp_mixed.reindex(labels=columns, axis='columns'), vmin=-12, vmax=np.max(Moods_means['p-value']),
            cmap='bwr_r', linewidths=0.35, linecolor='grey', center=np.log10(0.05), ax=axs[1,0], cbar=False)
sns.heatmap(bp_ADS.reindex(labels=columns, axis='columns'), vmin=-12, vmax=np.max(Moods_means['p-value']), cmap='bwr_r',
            linewidths=0.35, linecolor='grey', center=np.log10(0.05), ax=axs[1,1], cbar=False)
sns.heatmap(bp_min.reindex(labels=columns, axis='columns'), vmin=-12, vmax=np.max(Moods_means['p-value']), cmap='bwr_r',
            linewidths=0.35, linecolor='grey', center=np.log10(0.05), ax=axs[1,2], cbar=False)

sns.heatmap(prot_mixed.reindex(labels=columns, axis='columns'), vmin=-12, vmax=np.max(Moods_means['p-value']),
            cmap='bwr_r', linewidths=0.35, linecolor='grey', center=np.log10(0.05), ax=axs[2,0], cbar=False)
sns.heatmap(prot_ADS.reindex(labels=columns, axis='columns'), vmin=-12, vmax=np.max(Moods_means['p-value']),
            cmap='bwr_r', linewidths=0.35, linecolor='grey', center=np.log10(0.05), ax=axs[2,1], cbar=False)
sns.heatmap(prot_min.reindex(labels=columns, axis='columns'), vmin=-12, vmax=np.max(Moods_means['p-value']),
            cmap='bwr_r', linewidths=0.35, linecolor='grey', center=np.log10(0.05), ax=axs[2,2], cbar=False)

sns.heatmap(ratio_mixed.reindex(labels=columns, axis='columns'), vmin=-12, vmax=np.max(Moods_means['p-value']),
            cmap='bwr_r', linewidths=0.35, linecolor='grey', center=np.log10(0.05), ax=axs[3,0], cbar=False)
sns.heatmap(ratio_ADS.reindex(labels=columns, axis='columns'), vmin=-12, vmax=np.max(Moods_means['p-value']),
            cmap='bwr_r', linewidths=0.35, linecolor='grey', center=np.log10(0.05), ax=axs[3,1], cbar=False)
sns.heatmap(ratio_min.reindex(labels=columns, axis='columns'), vmin=-12, vmax=np.max(Moods_means['p-value']),
            cmap='bwr_r',linewidths=0.35, linecolor='grey', center=np.log10(0.05), ax=axs[3,2], cbar=False)

axes_all = [axs[0,0], axs[0,1], axs[0,2], axs[1,0], axs[1,1], axs[1,2], axs[2,0], axs[2,1], axs[2,2], axs[3,0],
            axs[3,1], axs[3,2]]

for ax in axes_all:
    ax.set_xticklabels(ax.get_xmajorticklabels(), size=32)
    ax.set_yticklabels(ax.get_ymajorticklabels(), size=28, rotation=45)
    ax.set_ylabel('')

axs[0,0].set_ylabel('Correlation of mutational effects', size=36, labelpad=8)
axs[1,0].set_ylabel('Correlation of mutational effects', size=36, labelpad=8)
axs[2,0].set_ylabel('Correlation of mutational effects', size=36, labelpad=8)
axs[3,0].set_ylabel('Correlation of mutational effects', size=36, labelpad=8)

axs[3,0].set_xlabel('Bias towards transcriptional mutations', size=36, labelpad=8)
axs[3,1].set_xlabel('Bias towards transcriptional mutations', size=36, labelpad=8)
axs[3,2].set_xlabel('Bias towards transcriptional mutations', size=36, labelpad=8)

axs[0,0].set_title('Mixed model', size=75, pad=75)
axs[0,1].set_title('ADS model', size=75, pad=75)
axs[0,2].set_title('Minimal model', size=75, pad=75)

cbar_ax.set_yticklabels(cbar_ax.get_yticklabels(), size=42)
cbar_ax.set_ylabel('Log10 p-value', size=54)

fig.subplots_adjust(left=0.195)

axs[0,0].annotate('Transcription\nrate', (-0.55, 0.5), xycoords='axes fraction', fontsize=65, annotation_clip=False,
                  horizontalalignment='center', verticalalignment='center')
axs[1,0].annotate('Translation\nrate', (-0.55, 0.5), xycoords='axes fraction', fontsize=65, annotation_clip=False,
                  horizontalalignment='center', verticalalignment='center')
axs[2,0].annotate('Protein\nabundance', (-0.55, 0.5), xycoords='axes fraction', fontsize=65, annotation_clip=False,
                  horizontalalignment='center', verticalalignment='center')
axs[3,0].annotate('Divergence\nratio', (-0.55, 0.5), xycoords='axes fraction', fontsize=65, annotation_clip=False,
                  horizontalalignment='center', verticalalignment='center')

fig.savefig('Moods_heatmap.pdf')
plt.close()

# 2) For KS statistics
KS_all = pd.DataFrame(columns=['Model', 'Comparison', 'Property', 'KS stats', 'KS p-value', 'Run', 'Mut_corr'])

for directory in iters:
    os.chdir(f'F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2/'
             f'Mut_correlations_WGD_1e5/{directory}')

    df_KS = pd.concat(map(pd.read_csv, glob.glob('KS_tests*.csv')))
    df_KS['Mut_corr'] = corr_vals[directory]

    KS_all = pd.concat([KS_all, df_KS])

os.chdir('F:/Simon_Aube/Hausser_2019_reanalysis/Simulations_ADS_Hausser/pythonProject/Final_simulations_v2/'
         'Mut_correlations_WGD_1e5')

# The 'No Cost' values are deleted
KS_all = KS_all[KS_all['Model'] != 'No Cost']

# 'iter{number}' are deleted from the Run tags
KS_all = KS_all.reset_index(drop=True)
for row in range(KS_all.shape[0]):
    KS_all.at[row, 'Run'] = KS_all.at[row, 'Run'][:-6]

# Mean KS statistics are computed
KS_bias = KS_all.groupby(by=['Model', 'Comparison', 'Property', 'Run', 'Mut_corr'], as_index=False).mean()

# And KS statistics are multiplied by -1
KS_bias['KS stats'] = KS_bias['KS stats'] * -1

KS_bias = KS_bias.astype({'Mut_corr':'float64'})

KS_comps = PdfPages(f'KS_bias_heatmap.pdf')

for comp in ['All duplicates', 'WGD', 'SSD']:
    comp_subset = KS_bias[KS_bias['Comparison'] == comp]

    # Subsets per parameter
    bm_KS = comp_subset[comp_subset['Property'] == 'Transcription rate'].copy()
    bp_KS = comp_subset[comp_subset['Property'] == 'Translation rate'].copy()
    prot_KS = comp_subset[comp_subset['Property'] == 'Protein abundance'].copy()
    div_KS = comp_subset[comp_subset['Property'] == 'Divergence ratio'].copy()

    # Subsets per model
    bm_mixed = bm_KS[bm_KS['Model'] == 'Mixed'].pivot(index='Mut_corr', columns='Run', values='KS stats')
    bp_mixed = bp_KS[bp_KS['Model'] == 'Mixed'].pivot(index='Mut_corr', columns='Run', values='KS stats')
    prot_mixed = prot_KS[prot_KS['Model'] == 'Mixed'].pivot(index='Mut_corr', columns='Run', values='KS stats')
    ratio_mixed = div_KS[div_KS['Model'] == 'Mixed'].pivot(index='Mut_corr', columns='Run', values='KS stats')

    bm_ADS = bm_KS[bm_KS['Model'] == 'ADS'].pivot(index='Mut_corr', columns='Run', values='KS stats')
    bp_ADS = bp_KS[bp_KS['Model'] == 'ADS'].pivot(index='Mut_corr', columns='Run', values='KS stats')
    prot_ADS = prot_KS[prot_KS['Model'] == 'ADS'].pivot(index='Mut_corr', columns='Run', values='KS stats')
    ratio_ADS = div_KS[div_KS['Model'] == 'ADS'].pivot(index='Mut_corr', columns='Run', values='KS stats')

    bm_min = bm_KS[bm_KS['Model'] == 'Minimal'].pivot(index='Mut_corr', columns='Run', values='KS stats')
    bp_min = bp_KS[bp_KS['Model'] == 'Minimal'].pivot(index='Mut_corr', columns='Run', values='KS stats')
    prot_min = prot_KS[prot_KS['Model'] == 'Minimal'].pivot(index='Mut_corr', columns='Run', values='KS stats')
    ratio_min = div_KS[div_KS['Model'] == 'Minimal'].pivot(index='Mut_corr', columns='Run', values='KS stats')

    # Column names are modified
    dfs = [bm_mixed, bp_mixed, prot_mixed, ratio_mixed, bm_ADS, bp_ADS, prot_ADS, ratio_ADS, bm_min, bp_min, prot_min,
           ratio_min]
    for df in dfs:
        df.columns = ['1/2', '1', '10', '2', '3', '4', '5', '6', '7', '8', '9']

    # For reordering the matrix
    columns = ['1/2', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

    # Generation of the figure
    fig, axs = plt.subplots(4, 3, figsize=[48, 48])

    fig.subplots_adjust(right=0.86)
    cbar_ax = fig.add_axes([0.9, 0.20, 0.02, 0.6])

    sns.heatmap(bm_mixed.reindex(labels=columns, axis='columns'), vmin=np.min(KS_bias['KS stats']),
                vmax=np.max(KS_bias['KS stats']), cmap='viridis',
                ax=axs[0, 0], cbar_ax=cbar_ax)
    sns.heatmap(bm_ADS.reindex(labels=columns, axis='columns'), vmin=np.min(KS_bias['KS stats']),
                vmax=np.max(KS_bias['KS stats']), cmap='viridis',
                ax=axs[0, 1], cbar=False)
    sns.heatmap(bm_min.reindex(labels=columns, axis='columns'), vmin=np.min(KS_bias['KS stats']),
                vmax=np.max(KS_bias['KS stats']), cmap='viridis',
                ax=axs[0, 2], cbar=False)

    sns.heatmap(bp_mixed.reindex(labels=columns, axis='columns'), vmin=np.min(KS_bias['KS stats']),
                vmax=np.max(KS_bias['KS stats']), cmap='viridis',
                ax=axs[1, 0], cbar=False)
    sns.heatmap(bp_ADS.reindex(labels=columns, axis='columns'), vmin=np.min(KS_bias['KS stats']),
                vmax=np.max(KS_bias['KS stats']), cmap='viridis',
                ax=axs[1, 1], cbar=False)
    sns.heatmap(bp_min.reindex(labels=columns, axis='columns'), vmin=np.min(KS_bias['KS stats']),
                vmax=np.max(KS_bias['KS stats']), cmap='viridis',
                ax=axs[1, 2], cbar=False)

    sns.heatmap(prot_mixed.reindex(labels=columns, axis='columns'), vmin=np.min(KS_bias['KS stats']),
                vmax=np.max(KS_bias['KS stats']), cmap='viridis',
                ax=axs[2, 0], cbar=False)
    sns.heatmap(prot_ADS.reindex(labels=columns, axis='columns'), vmin=np.min(KS_bias['KS stats']),
                vmax=np.max(KS_bias['KS stats']), cmap='viridis',
                ax=axs[2, 1], cbar=False)
    sns.heatmap(prot_min.reindex(labels=columns, axis='columns'), vmin=np.min(KS_bias['KS stats']),
                vmax=np.max(KS_bias['KS stats']), cmap='viridis',
                ax=axs[2, 2], cbar=False)

    sns.heatmap(ratio_mixed.reindex(labels=columns, axis='columns'), vmin=np.min(KS_bias['KS stats']),
                vmax=np.max(KS_bias['KS stats']), cmap='viridis',
                ax=axs[3, 0], cbar=False)
    sns.heatmap(ratio_ADS.reindex(labels=columns, axis='columns'), vmin=np.min(KS_bias['KS stats']),
                vmax=np.max(KS_bias['KS stats']), cmap='viridis',
                ax=axs[3, 1], cbar=False)
    sns.heatmap(ratio_min.reindex(labels=columns, axis='columns'), vmin=np.min(KS_bias['KS stats']),
                vmax=np.max(KS_bias['KS stats']), cmap='viridis',
                ax=axs[3, 2], cbar=False)

    axes_all = [axs[0, 0], axs[0, 1], axs[0, 2], axs[1, 0], axs[1, 1], axs[1, 2], axs[2, 0], axs[2, 1], axs[2, 2],
                axs[3, 0],
                axs[3, 1], axs[3, 2]]

    for ax in axes_all:
        ax.set_xticklabels(ax.get_xmajorticklabels(), size=32)
        ax.set_yticklabels(ax.get_ymajorticklabels(), size=28, rotation=45)
        ax.set_ylabel('')

    axs[0, 0].set_ylabel('Correlation of mutational effects', size=36, labelpad=8)
    axs[1, 0].set_ylabel('Correlation of mutational effects', size=36, labelpad=8)
    axs[2, 0].set_ylabel('Correlation of mutational effects', size=36, labelpad=8)
    axs[3, 0].set_ylabel('Correlation of mutational effects', size=36, labelpad=8)

    axs[3, 0].set_xlabel('Bias towards transcriptional mutations', size=36, labelpad=8)
    axs[3, 1].set_xlabel('Bias towards transcriptional mutations', size=36, labelpad=8)
    axs[3, 2].set_xlabel('Bias towards transcriptional mutations', size=36, labelpad=8)

    axs[0, 0].set_title('Mixed model', size=75, pad=75)
    axs[0, 1].set_title('ADS model', size=75, pad=75)
    axs[0, 2].set_title('Minimal model', size=75, pad=75)

    cbar_ax.set_yticklabels(cbar_ax.get_yticklabels(), size=42)
    cbar_ax.set_ylabel('-(KS statistics)', size=54, labelpad=10)

    fig.subplots_adjust(left=0.195)

    axs[0, 0].annotate('Transcription\nrate', (-0.55, 0.5), xycoords='axes fraction', fontsize=65,
                       annotation_clip=False,
                       horizontalalignment='center', verticalalignment='center')
    axs[1, 0].annotate('Translation\nrate', (-0.55, 0.5), xycoords='axes fraction', fontsize=65, annotation_clip=False,
                       horizontalalignment='center', verticalalignment='center')
    axs[2, 0].annotate('Protein\nabundance', (-0.55, 0.5), xycoords='axes fraction', fontsize=65, annotation_clip=False,
                       horizontalalignment='center', verticalalignment='center')
    axs[3, 0].annotate('Divergence\nratio', (-0.55, 0.5), xycoords='axes fraction', fontsize=65, annotation_clip=False,
                       horizontalalignment='center', verticalalignment='center')

    if comp == 'All duplicates':
        comp = 'All duplicate'

    fig.suptitle(f'Comparisons with the empirical distributions for the three models\n({comp} couples considered)',
                 size=85)

    fig.savefig(KS_comps, format='pdf')

    fig.clf()

KS_comps.close()

# 3) Mean KS statistics
mean_KS = KS_all.groupby(by=['Model', 'Comparison', 'Run', 'Mut_corr'], as_index=False).mean()
mean_KS['KS stats'] = mean_KS['KS stats'] * -1

# Only comparison with all duplicate couples are kept
mean_KS = mean_KS[mean_KS['Comparison'] == 'All duplicates']
mean_KS = mean_KS.astype({'Mut_corr':'float64'})

# Matrices for each model
KS_mixed = mean_KS[mean_KS['Model'] == 'Mixed'].pivot(index='Mut_corr', columns='Run', values='KS stats')
KS_ADS = mean_KS[mean_KS['Model'] == 'ADS'].pivot(index='Mut_corr', columns='Run', values='KS stats')
KS_min = mean_KS[mean_KS['Model'] == 'Minimal'].pivot(index='Mut_corr', columns='Run', values='KS stats')

# Column names are modified
dfs = [KS_mixed, KS_ADS, KS_min]
for df in dfs:
    df.columns = ['1/2', '1', '10', '2', '3', '4', '5', '6', '7', '8', '9']

# For reordering the matrix
columns = ['1/2', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

# The figure is constructed
fig, axs = plt.subplots(1, 3, figsize=[45, 15])

fig.subplots_adjust(right=0.86)
cbar_ax = fig.add_axes([0.9, 0.20, 0.02, 0.6])

sns.heatmap(KS_mixed.reindex(labels=columns, axis='columns'), vmin=np.min(mean_KS['KS stats']),
            vmax=np.max(mean_KS['KS stats']), cmap='viridis', square=True,
            ax=axs[0], cbar_ax=cbar_ax)
sns.heatmap(KS_ADS.reindex(labels=columns, axis='columns'), vmin=np.min(mean_KS['KS stats']),
            vmax=np.max(mean_KS['KS stats']), cmap='viridis', square=True,
            ax=axs[1], cbar=False)
sns.heatmap(KS_min.reindex(labels=columns, axis='columns'), vmin=np.min(mean_KS['KS stats']),
            vmax=np.max(comp_subset['KS stats']), cmap='viridis', square=True,
            ax=axs[2], cbar=False)

for ax in [axs[0], axs[1], axs[2]]:
    ax.set_xticklabels(ax.get_xmajorticklabels(), size=32)
    ax.set_yticklabels(ax.get_ymajorticklabels(), size=28, rotation=45)
    ax.set_ylabel('')

axs[0].set_ylabel('Correlation of mutational effects', size=36, labelpad=8)

axs[0].set_xlabel('Bias towards transcriptional mutations', size=36, labelpad=8)
axs[1].set_xlabel('Bias towards transcriptional mutations', size=36, labelpad=8)
axs[2].set_xlabel('Bias towards transcriptional mutations', size=36, labelpad=8)

axs[0].set_title('Mixed model', size=75, pad=55)
axs[1].set_title('ADS model', size=75, pad=55)
axs[2].set_title('Minimal model', size=75, pad=55)

cbar_ax.set_yticklabels(cbar_ax.get_yticklabels(), size=42)
cbar_ax.set_ylabel('-(KS statistics)', size=54, labelpad=10)

fig.savefig('Mean_Ks_heatmap.pdf')
plt.close()
