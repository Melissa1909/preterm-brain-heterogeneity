# Helpers for plotting analysis results

import os
from os.path import join
import pandas as pd
import numpy as np
from enigmatoolbox.utils.parcellation import parcel_to_surface
from enigmatoolbox.plotting import plot_cortical
from scipy.stats import spearmanr

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter,FuncFormatter
import seaborn as sns
from statannotations.Annotator import Annotator
import sys
sys.path.append('')
from spatial_heterogeneity import calc_infra_supra_percentage


def get_significant_cortical_vals_for_plotting(data, rois_cortical, statistic='t_statistic', bilateral=True):
    '''
    Prepare data for plotting significant cortical values.

    data: pandas dataframe with t-statistics and p-values for each ROI
    rois_cortical: list of ROIs to be plotted
    statistic: statistic to be plotted, default='t_statistic', has to be in "data"
    bilateral: if True, duplicate values to receive 68 values necessary for plotting with ENIGMA toolbox
    '''
    # get significant rois
    sig = data.loc[data['ROI'].isin(rois_cortical)].copy()
    sig.loc[sig['p_fdr'] > 0.05, sig.columns[1]] = 0  # set non-significant rois' t-statistics to 0
    
    sig_vals = sig[[statistic]].to_numpy().squeeze()
    
    if bilateral == True:
        # duplicate values to receive 68 values necessary for plotting with ENIGMA toolbox
        sig_vals = np.tile(sig_vals,(1,2)).squeeze()
    
    return sig_vals


def plot_brain_map(data, outname, cmap='RdBu_r', limits=(-5,5), scale=1, fill=0):
    '''
    Plot a cortical rendering of the values in data.

    data: numpy array of values per ROI in Desikan-Killiany atlas, needs to be 68 values
    outname: file where figure should be saved to
    scale: image quality, default=1
    limit: color range limit
    '''
    # # if data is not 68 values, raise an error
    # if len(data) != 68:
    #     raise ValueError('Data needs to be 68 values for plotting with ENIGMA toolbox.')
    
    # convert data to surface
    surf = parcel_to_surface(data, 'aparc_fsa5', fill=fill)

    # plot and save
    plot_cortical(array_name=surf, 
                    surface_name="fsa5", 
                    size=(1000, 300),
                    cmap=cmap, 
                    color_bar=True, 
                    color_range=(limits[0], limits[1]),
                    screenshot=True, 
                    filename=outname,
                    scale=(scale, scale), 
                    dpi=300)
    print(f'Plotted brain map and saved it to {outname}')


def plot_group_difference_global(data, x='dx', y='WMV', yaxis_label = f'WMV [mm$^3$]', 
                                    outname='group_difference_WMV.svg'):
    '''
    Plot violin plot of global measures (e.g. WMV, GMV, mean CT) for term and preterm groups.

    data: pd.DataFrame, df containing the data for each subject
    x: str, group variable (e.g. 'dx')
    y: str, variable to plot (e.g. 'WMV')
    yaxis_label: str, label for the y-axis
    outname: str, output filename for the figure
    '''
    #cm = 1/2.54  # cm to inch
    order = ['CN', 'preterm']
    plt.figure(figsize=(3,4))
    ax = sns.violinplot(data=data, x=x, y=y, order=order, palette=['royalblue', 'darkorange'], 
                        fill=True, inner='box', linewidth=0.5)
    ax.set_xlabel('')
    ax.set_xticklabels(['Term', 'Preterm'])
    ax.set_ylabel(yaxis_label)
    
    # format y axis
    # ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x/1e5:.1f}e5'))
    ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    sns.despine()

    # add significance annotations
    pairs=[('CN', 'preterm')]
    annotator = Annotator(ax, pairs, data=data, x=x, y=y, order=order)
    annotator.configure(test='t-test_ind', text_format='star', loc='inside', comparisons_correction='fdr_bh', line_width=0.5)
    annotator.apply_and_annotate()

    # save figure
    plt.savefig(outname, bbox_inches='tight')
    plt.show()
    print(f'Plotted group difference for {y} and saved it to {outname}')


def prepare_data_for_brain_plot(cortical_data, sub_id, filter_attribute, bilateral=True):
    '''
    Duplicate and slice data for plotting brain maps with plot_brain_map.

    cortical_data: pd.DataFrame that will be filtered for filter_attribute
    sub_id: str, subject id
    filter_attribute: str, attribute to filter data for (e.g., centile_CT)
    bilateral: bool, if True, duplicate data to receive 68 values, default=True
    '''
    cortical_data_sub = cortical_data[cortical_data['participant']==sub_id]  # only extract subject's data
    cortical_data_sub = cortical_data_sub.filter(regex=filter_attribute)  # only extract the attribute of interest
    cortical_data_sub = cortical_data_sub.to_numpy().flatten().astype(float)
    if bilateral == True:  # duplicate data if necessary to receive 68 values
        cortical_data_sub = np.tile(cortical_data_sub, (1,2)).squeeze()
    return cortical_data_sub


def plot_individual_brain_maps(cortical_data, sub_id, outdir, scale=1, cmap='YlOrRd', limits=(0, 1), 
                                infra_supra=False, brain_measure='CT', bilateral=True):
    '''
    Function to plot individual subjetct's deviation maps.

    cortical_data: pd.df, pandas dataframe with rCTD values
    sub_id: str, subject id, e.g., BEST-BN-001
    outdir: os.path, directory where figures should be saved to
    scale: int, scale of the figure, default=1
    cmap: colormap, default='YlOrRd'
    limits: tuple, color range limits, default=(0, 1)
    infra_supra: bool, if True, plot infra- and supranormal values only in a separate figure
    brain_measure: str, measure to be plotted, default='CT'
    bilateral: bool, if True, duplicate data to receive 68 values, default=True
    '''
    # convert rCTD to numpy and filter necessary data
    cortical_data_sub = prepare_data_for_brain_plot(cortical_data, sub_id, filter_attribute=f'centile_{brain_measure}')

    plot_brain_map(cortical_data_sub, outname=os.path.join(outdir, f'{sub_id}_{brain_measure}_individual_centiles.svg'), 
                    cmap=cmap, limits=limits, scale=scale, fill=0.5)
    print(f'Plotted individual rCTD map for subject {sub_id} and saved it to {outdir}')
    
    if infra_supra == True:
        # plot infra/supra values only
        sig_deviations = np.where(cortical_data_sub < 0.05, -1, np.where(cortical_data_sub > 0.95, 1, 0))
        outname_is = os.path.join(outdir, f'{sub_id}_{brain_measure}_individual_centiles_extranormal.svg')

        plot_brain_map(sig_deviations, outname=outname_is, 
                    cmap='RdBu_r', limits=(-1,1), scale=scale, fill=0)
    

def correlation_plot(x, y, data, color, xlabel, ylabel, outname, one_sided=False):
    '''
    Function to plot relation of x and y.
    
    x: str, variable name in data
    y: str, variable name in data
    data: pd.Dataframe, df with variables x and y
    color: str, color of the plot
    xlabel: str, x-axis label
    ylabel: str, y-axis label
    outname: os.path, name of the file where figure should be saved
    '''
    # Spearman correlation between x and y
    r, p = spearmanr(data[x], data[y], nan_policy='omit')
    if one_sided == True:
        p = p/2

    # plotting
    matplotlib.rcParams['grid.linewidth'] = 6
    cm = 1/2.54
    fig = plt.figure(figsize=(5.5*cm,4*cm), dpi=150)
    
    plot = sns.regplot(x=data[x], y=data[y], color=color, n_boot=10000, 
                        scatter_kws={'s': 3, 'edgecolor': 'None'}, truncate=True, line_kws={'linewidth': 0.75})
    
    # add correlation coefficient and p-value
    textstr = '\n'.join((
            r'$r=%.3f$' % (r),
            r'$p=%.3f$' % (p)))
    if one_sided == True:
            textstr = '\n'.join((
                r'$r=%.3f$' % (r),
                r'$p_1=%.3f$' % (p)))
    props = dict(boxstyle='round', facecolor='silver', alpha=0.4, edgecolor='silver')
    plot.text(0.03, 0.98, textstr, transform=plot.transAxes, fontsize=5, verticalalignment='top', bbox=props)

    # label settings
    plt.xlabel(xlabel, fontsize=6)
    plt.xticks(fontsize=5)
    plt.ylabel(ylabel, fontsize=6)
    plt.yticks(fontsize=5)
    plt.tick_params(axis='both', which='major', width=0.25, pad=1)
    
    # add xticklabels for SES
    if x == "SES_at_birth":
        plt.xticks([1,2,3],['high', 'middle', 'low'])
    
    plt.tight_layout()
    sns.despine()
    plt.savefig(outname, dpi=300)
    plt.show()
    
    
def plot_percentage_outside_norm_groups(rois, ft_scores_df, pt_scores_df, age_thr_lower, age_thr_upper, outname, add_legend=False):
    '''
    Plot the percentage outside the normal range (i.e., <5th or >95th percentile) for each ROI.

    rois: list of ROIs
    ft_scores_df: pd.DataFrame, deviation scores of FT group
    pt_scores_df: pd.DataFrame, deviation scores of PT group
    age_thr_lower: str, lower age threshold of gestational age
    age_thr_upper: str, upper age threshold of gestational age
    outname: os.path, name of the file where figure should be saved
    add_legend: bool, whether to add a legend, default=False
    '''
    rois_cent = ['centile_' + roi for roi in rois]
    infra_supra_all = calc_infra_supra_percentage(rois_cent, ft_scores_df, pt_scores_df)

    fsize=8
    ax = infra_supra_all.plot(x='rois', y=['infranormal_ft','supranormal_ft','infranormal_pt','supranormal_pt'], kind="barh", rot=0,
                                figsize=(3,12),fontsize=fsize,color=['dimgray','lightgray', 'mediumblue', 'red'],
                                legend=False, width=0.75)
    # label settings
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('Subjects with extranormal deviation [%]',fontsize=fsize+2)
    ax.set_ylabel('')
    if rois[0].startswith('CT_'):
        rois_cortical_names = ['CTh ' + roi.split('_')[1] for roi in rois]
    else:
        rois_cortical_names = ['SA ' + roi.split('_')[1] for roi in rois]
    ax.set_yticklabels(rois_cortical_names)
        
    # set the x-axis limit    
    plt.xlim(0, 40)

    # add legend
    if add_legend is True:
        legend=ax.legend(['Infranormal FT','Supranormal FT','Infranormal PT','Supranormal PT'], fontsize=fsize, loc ='lower right')
        legend.get_frame().set_facecolor('white')

    sns.despine()
    
    plt.title('GA between {0:.0f} and {1:.0f}'.format(age_thr_lower, age_thr_upper), fontsize=fsize+2)
    plt.savefig(outname, bbox_inches="tight", dpi=600)
    plt.show()
    
    return infra_supra_all


def plot_celltype_correlations(data, outname, xlabels_bool=False):
    cm = 1/2.54
    fig = plt.figure(figsize=(14 * cm, 4.5 * cm), dpi=300)
    
    ax = sns.heatmap(data, cmap='coolwarm', xticklabels=xlabels_bool, yticklabels=False, 
                    cbar_kws={'label': 'Spearman rho', 'shrink': 0.8})
    
    if xlabels_bool == True:
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=5, rotation=0)
        plt.tick_params(axis='x', length=0)
    
    # Customize the colorbar's font size
    cbar = ax.collections[0].colorbar  # Get the colorbar
    cbar.ax.tick_params(labelsize=5)  # Set font size for tick labels
    cbar.ax.set_ylabel('Spearman rho', fontsize=5)  # Set font size for the label
    cbar.ax.tick_params(width=0.25) 
    
    #plt.tight_layout()
    plt.savefig(outname, dpi=600)
    plt.show()
    
    

def plot_gene_exp(celltype, it_dir, out_dir, scale=1):
    '''
    Function to plot gene expression maps.

    celltype: celltype, e.g., 'Astro'
    it_dir: directory where input data is stored (i.e., expression values for that celltype calculated before)
    out_dir: directory where figure should be saved to
    scale: scale of the figure, default=1
    '''
    ge = pd.read_csv(join(it_dir, f'expression_{celltype}.csv'), delimiter=',', index_col=0)
    ge = ge[celltype].values.flatten().astype(float)
    ### BUG in previous version; used np.tile
    
    outname = os.path.join(out_dir, f'{celltype}_expression.svg')
    plot_brain_map(ge, outname, cmap='YlOrRd', limits=(0.4, 0.7), scale=scale, fill=0.5)








# def plot_percentage_outside_norm_supp(rois, ft_scores_df, pt_scores_df, age_thr_lower, age_thr_upper, outname, add_legend=False):
#     '''
#     Plot the percentage outside the normal range (i.e., <5th or >95th percentile) for each ROI.

#     rois: list of ROIs
#     ft_scores_df: dev scores of FT group
#     pt_scores_df: dev scores of PT group
#     age_thr_lower: what was used as a age threshold; e.g., 28 weeks
#     age_thr_upper: what was used as a age threshold; e.g., 28 weeks
#     outname: name of the file where figure should be saved
#     add_legend: whether to add a legend, default is False
#     '''
#     infra_supra_all = calc_infra_supra_amounts(rois, ft_scores_df, pt_scores_df)

#     fsize=8
#     ax = infra_supra_all.plot(x='rois', y=['infranormal_ft','supranormal_ft','infranormal_pt','supranormal_pt'], kind="barh", rot=0,
#                                 figsize=(3,12),fontsize=fsize,color=['dimgray','lightgray', 'mediumblue', 'red'],
#                                 legend=False)

#     ax.invert_yaxis()  # labels read top-to-bottom
#     ax.set_xlabel('Subjects with extranormal deviation [%]',fontsize=8)
#     ax.set_ylabel('')
#     rois_cortical_names = ['CTh ' + roi.split('_')[1] for roi in rois]
#     ax.set_yticklabels(rois_cortical_names)

#     # set the x-axis limit
#     #xlim = np.round(infra_supra_all[['infranormal_ft','supranormal_ft','infranormal_pt','supranormal_pt']].max().max()*1.1)
#     #print(xlim)
#     xlim = 60
#     plt.xlim(0, xlim)

#     # add legend
#     if add_legend is True:
#         legend=ax.legend(['Infranormal FT','Supranormal FT','Infranormal PT','Supranormal PT'],fontsize=fsize, loc ='lower right')
#         legend.get_frame().set_facecolor('white')


#     plt.title('GA between {0} and {1}'.format(age_thr_lower, age_thr_upper), fontsize=fsize+2)
#     plt.savefig(outname,bbox_inches="tight",dpi=300)
#     plt.show()
