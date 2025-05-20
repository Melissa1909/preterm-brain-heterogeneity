# Analysis helpers for braincharts processing

import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
sns.set_style('white')
from statsmodels.stats.multitest import multipletests
from statsmodels.formula.api import ols
import math


def calculate_cohens_d(data, roi):
    '''
    Calculates Cohen's d to estimate effect sizes.
    data: pandas df
    roi: name of roi, e.g. 'lh_CT_bankssts'
    '''
    
    g1 = data[data['diagnosis'] == 0]
    g2 = data[data['diagnosis'] == 1]
    
    # preterm stats
    mean1 = g1[roi].mean()
    std1 = g1[roi].std()
    n1 = g1[roi].count()
        
    # fullterm stats
    mean2 = g2[roi].mean()
    std2 = g2[roi].std()
    n2 = g2[roi].count()
        
    # cohen's d calculation
    pooled_std = math.sqrt(((n1-1)*std1**2 + (n2-1)*std2**2) / (n1 + n2 - 2))
    cohen_d = (mean1 - mean2) / pooled_std
    
    return cohen_d


def group_comparison(rois, dat, covariates=['sex', 'age_days']):
    '''
    Perform a group comparison for each roi in rois using a linear model with the diagnosis as the predictor and the covariates as confounders.
    Will also perform multiple comparisons correction using the Benjamini-Hochberg method.

    rois: list of ROIs to perform the group comparison
    dat: DataFrame containing the data
    covariates: list of covariates to include in the model
    '''
    # Initialize an empty list to store the results
    results = []

    for roi in rois:
        # add binary variable for diagnosis
        dat['diagnosis'] = dat['dx'].map({'preterm': 1, 'CN': 0})

        # fit the model
        formula = f'{roi} ~ diagnosis + {" + ".join(covariates)}'
        model = ols(formula, data=dat).fit()
        
        # extract the t-value and p-value for the diagnosis variable
        t_value = model.tvalues['diagnosis']
        p_value = model.pvalues['diagnosis']
        
        # calculate Cohen's d
        d = calculate_cohens_d(dat, roi)
        
        results.append((roi, t_value, d, p_value))

    # convert the results list to a DataFrame
    result_df = pd.DataFrame(results, columns=['ROI', 't_statistic', 'Cohen_d', 'p_value'])

    # correct for multiple comparisons
    _, p_fdr, _, _ = multipletests(result_df['p_value'], method='fdr_bh')
    result_df['p_fdr'] = p_fdr

    return result_df


def get_centile_scores_per_subject(analysis_dir, rois, base_file='result_deviation_CT_bankssts.csv'):
    '''
    Load centile scores for each ROI computed in the braincharts framework.
    
    analysis_dir: os.path, path to the directory where the braincharts results are stored
    rois: list of str, list of ROIs for which centile scores are computed
    base_file: str, name of one of the braincharts results files to load the other necessary information
    '''

    # Load the base dataframe
    base_df = pd.read_csv(os.path.join(analysis_dir, base_file))

    # Extract participant IDs
    participants = base_df['participant'].values
    session = base_df['session'].values

    # Initialize an empty array to store centile scores
    centile_scores = np.zeros((len(participants), len(rois)))

    # # Iterate over each ROI to collect centile scores
    for i, roi in enumerate(rois):
        # Skip entorhinal modeling for SA
        if roi == 'SA_entorhinal':
            centile_scores[:, i] = np.nan
            continue
        
        # Load the centile score data for the current ROI
        roi_file = os.path.join(analysis_dir, f'result_deviation_{roi}.csv')
        roi_data = pd.read_csv(roi_file)
        
        # Align the centile scores with the participant IDs
        centile_scores[:, i] = roi_data['centile_score'].values

    # Convert the array to a DataFrame
    cent_df = pd.DataFrame(centile_scores, columns=[f'centile_{roi}' for roi in rois])

    # Add the participant IDs back to the DataFrame
    cent_df.insert(0, 'participant', participants)
    cent_df.insert(1, 'session', session)

    return cent_df


def get_amount_infra_supra(cent_scores, infra_thr = 0.05, supra_thr = 0.95):
    '''
    Count infranormal and supranormal per subject.

    cent_scores: dataframe with centile scores for each ROI
    infra_thr: float, threshold for infranormal values, default: 0.05
    supra_thr: float, threshold for supranormal values, default: 0.95
    '''
    cent_scores['amount_infranormal'] = (cent_scores < infra_thr).sum(axis=1)
    cent_scores['amount_supranormal'] = (cent_scores > supra_thr).sum(axis=1)
    cent_scores.reset_index(inplace=True)
    return cent_scores


def percent_extranormal(df):
    '''
    Calculate the percentage of subjects with extranormal values.
    '''
    return df[(df['amount_infranormal']>0) | (df['amount_supranormal']>0)].shape[0] / df.shape[0] * 100


def calc_infra_supra_percentage(rois, ft_scores_df, pt_scores_df, infra_thr = 0.05, supra_thr = 0.95):
    '''
    Calculate the amount of infra- and supranormal values for each ROI.

    rois: list of str, ROIs of deviation scores (e.g., centile_)
    ft_scores_df: pd.DataFrame, deviation scores of FT group
    pt_scores_df: pd.DataFrame, deviation scores of PT group
    '''
    # initialize result arrays
    infranormal_pt = np.zeros([len(rois),1])
    infranormal_ft = np.zeros([len(rois),1])
    supranormal_pt = np.zeros([len(rois),1])
    supranormal_ft = np.zeros([len(rois),1])

    for r, roi in enumerate(rois):

        # preterm: count how many values are <5th or >95th percentile
        data = pt_scores_df[roi]
        infranormal_pt[r,] = (data[data < infra_thr].count() ) / len(data) *100
        supranormal_pt[r,] = (data[data > supra_thr].count() ) / len(data) *100

        # fullterm
        data = ft_scores_df[roi]
        infranormal_ft[r,] = ( data[data < infra_thr].count() ) / len(data) *100
        supranormal_ft[r,] = ( data[data > supra_thr].count() ) / len(data) *100

    # transform to df
    infra_supra_all = pd.DataFrame(infranormal_pt, index=rois, columns=['infranormal_pt'])
    infra_supra_all['supranormal_pt'] = supranormal_pt
    infra_supra_all['infranormal_ft'] = infranormal_ft
    infra_supra_all['supranormal_ft'] = supranormal_ft
    infra_supra_all.reset_index(inplace=True)
    infra_supra_all.rename(columns={'index':'rois'},inplace=True)

    return infra_supra_all


def binarize_extranormal(data):
    '''
    Binarize extranormal deviations so that infra- or supranormal values are coded as 1, otherwise 0.
    
    data: df with deviation scores
    '''
    data_bin = pd.DataFrame(np.where((data < 0.05) | (data > 0.95), 1, 0), columns=data.columns)
    return data_bin


def plot_binarized_extranormal(data_bin, title, outname, bin=True):
    '''
    Plot the binarized extranormal deviations.
    
    data_bin: pd.DataFrame, binarized deviation scores of preterm subjects
    title: str, title of the plot
    outname: str, path where the plot should be saved to
    bin: bool, if True, use binary color map, otherwise use diverging color map
    '''
    if bin:
        fig = plt.figure(figsize=(2.5, 4))
        heatmap = sns.heatmap(data_bin, cmap=['white','k'], cbar=False, figure=fig)
        plt.xlabel('ROIs', fontsize=10)
        plt.ylabel('Subject-IDs', fontsize=10)
        
    else:  # for supplementary figure
        fig = plt.figure(figsize=(5, 2.5))
        heatmap = sns.heatmap(data_bin.T, cmap='RdBu_r', center=0.5, cbar=True, 
                                cbar_kws={'label': 'Deviation score'}, figure=fig)
        heatmap.collections[0].colorbar.ax.yaxis.label.set_size(10)   
        plt.xlabel('Subject-IDs', fontsize=10)
        plt.ylabel('ROIs', fontsize=10)     
        
    plt.title(title)
    heatmap.set_xticklabels('')
    heatmap.set_yticklabels('')

    plt.tight_layout()
    plt.savefig(outname, dpi=300)
    plt.show()


def plot_correlation_matrix_kde(corr_matrix, outname, fsize=12):
    '''
    Plot a heatmap of the correlation matrix and a KDE plot of the average Spearman rho across subjects.
    
    corr_matrix: pd.DataFrame, correlation matrix subject x subject
    outname: str, file where figure should be saved to
    fsize: int, fontsize
    '''
    # create figure with a grid for heatmap, kde, and color bar
    fig = plt.figure(figsize=(8, 6))
    gs = GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[15, 0.5], figure=fig)
    ax_heatmap = fig.add_subplot(gs[0, 0])
    ax_cbar = fig.add_subplot(gs[1, 0])
    
    # heatmap of cross-correlation between subjects
    heatmap = sns.heatmap(corr_matrix, cmap='vlag', vmax=1, vmin=-1, square=True, 
                            ax=ax_heatmap, cbar_ax=ax_cbar, 
                            cbar_kws={'orientation': 'horizontal', 'label': 'Spearman rho'}) 
    ax_heatmap.set_xlabel('Subject-IDs', fontsize=fsize)
    ax_heatmap.set_ylabel('Subject-IDs', fontsize=fsize)
    ax_heatmap.set_xticklabels('')
    ax_heatmap.set_yticklabels('')
    
    # KDE plot of avg spearman between subjects
    sns.set_style('whitegrid')
    ax_kde = fig.add_subplot(gs[0, 1])
    row_means = np.mean(corr_matrix, axis=1)  # average correlation per subject with all others
    sns.kdeplot(row_means, ax=ax_kde, vertical=True, fill=False, color='k')

    ax_kde.set_ylabel('Average Spearman rho across subjects', fontsize=fsize)
    ax_kde.set_xlabel('Density', fontsize=fsize)
    ax_kde.set_ylim(-1, 1)
    
    plt.tight_layout()
    plt.savefig(outname, dpi=300)
    plt.show()
    
    
def plot_mean_rho(corr_matrix, title, outname):
    '''
    Plot the distribution of mean Spearman rho.
    
    corr_matrix: pd.DataFrame, correlation matrix subject x subject
    title: str, title of the plot
    outname: str, path where the plot should be saved to
    '''
    row_means = corr_matrix.mean(axis=1)
    
    fsize=12
    sns.set_style('whitegrid')
    plt.figure(figsize=(5,2))
    
    sns.kdeplot(x=row_means, linewidth=2, color='k')
    plt.xlim(-1,1)
    sns.despine()
    plt.title(title, fontsize=fsize)
    plt.xlabel('Spearman rho', fontsize=fsize)
    plt.ylabel('Density', fontsize=fsize)
    plt.xticks(fontsize=fsize-2)
    plt.yticks(fontsize=fsize-2)

    plt.tight_layout()
    plt.savefig(outname, dpi=300)
    plt.show()