import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import pingouin
from statsmodels.stats.multitest import multipletests

import sys
sys.path.append('../code/')
from plotting import plot_individual_brain_maps, prepare_data_for_brain_plot, plot_brain_map


def plot_longitudinal_data_subplots(df, rois, modality, out_dir_fig=os.getcwd(), filter_extranormal=False, fontsize=18, num_plots=34, num_cols=7, dataset='BLS'):
    '''
    Plot longitudinal data for global measures in the rois list in individual subplots.
    
    df: DataFrame containing data for two timepoints
    rois: List of ROIs to plot (should be 34 in this case)
    modality: Modality of the data (e.g., CT)
    out_dir_fig: Output directory to save the figure
    filter_extranormal: Boolean flag to filter and plot only subjects with infra- or supranormal values
    '''
    # Define plot layout
    num_rows = (num_plots + num_cols - 1) // num_cols 

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(25/7*num_cols, num_rows*5))
    axes = axes.flatten()  # flatten to easily iterate over it


    for roi, (roi, ax) in enumerate(zip(rois, axes)):
        # extract roi name base
        roi_name = roi.split('_')[-1]

        # skip SA entorhinal centile plot
        if roi == 'centile_SA_entorhinal':
            ax.set_title(roi_name, fontsize=fontsize)
            continue

        # fill background for extranormal regions
        ax.axhspan(0, 0.05, facecolor='royalblue', alpha=0.4) 
        ax.axhspan(0.95, 1, facecolor='red', alpha=0.4)  
        
        for sub in df['participant'].unique():
            sub_data = df[df['participant'] == sub]

            # if filtering for extranormal values
            if filter_extranormal:
                infranormal = (sub_data[roi] < 0.05).any()
                supranormal = (sub_data[roi] > 0.95).any()
                
                if not (infranormal or supranormal):
                    continue

            # scatter plot points for each timepoint
            ax.scatter(x=sub_data['timepoint'], y=sub_data[roi], color='k', s=40, alpha=.8)
            # Draw a line connecting the timepoints for the same person
            ax.plot(sub_data['timepoint'], sub_data[roi], alpha=.8, color='k', linewidth=1.5)
            
            ax.set_title(roi_name, fontsize=fontsize)
            ax.set_ylim(-0.05, 1.05)
            ax.grid(False)
            ax.set_xticks([1, 2])
            ax.tick_params(axis='both', which='major', labelsize=fontsize-2)

    # If there are any unused subplots, remove them (only for non-global measures to keep size consistent)
    # if modality!='global':
    for i in range(num_plots, len(axes)):
        fig.delaxes(axes[i])

    # Choose filename based on whether filtering is applied
    if filter_extranormal:
        filename = f'{dataset}_{modality}_combined_longitudinal_extranormal_plotnum{num_plots}.svg'
    else:
        filename = f'{dataset}_{modality}_combined_longitudinal.svg'
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir_fig, filename))
    plt.show()



def save_source_data_for_extranormal(df, rois, outname):
    """
    Save source data for extranormal participants.
    
    df: DataFrame containing data for two timepoints
    rois: List of ROIs to include in the source data
    outname: Output filename for the source data CSV
    """
    # Generate a new df
    df_copy = df.copy()
    for roi in rois:
        df_copy[roi] = np.nan

    # keep track of included participants (i.e. those with at least one extranormal value for that ROI)
    included_participants = set()

    for roi in rois:
        if roi == 'centile_SA_entorhinal':  # skip as per plotting code
            continue

        # get mask of extranormal participants for this ROI
        is_extranormal = df.groupby('participant')[roi].apply(lambda x: (x < 0.05).any() or (x > 0.95).any())
        extranormal_ids = is_extranormal[is_extranormal].index.tolist()

        # mark these participants as included
        included_participants.update(extranormal_ids)

        # fill values only for extranormal participants
        for pid in extranormal_ids:
            df_copy.loc[df_copy['participant'] == pid, roi] = df.loc[df['participant'] == pid, roi]

    # filter to only included participants
    df_filtered = df_copy[df_copy['participant'].isin(included_participants)].copy()

    # #  anonymize participant IDs
    # unique_ids = df_filtered['participant'].unique()
    # id_map = {pid: f"sub-{i:03d}" for i, pid in enumerate(unique_ids, 1)}
    # df_filtered['participant'] = df_filtered['participant'].map(id_map)

    # save
    df_filtered = df_filtered[['participant', 'timepoint'] + rois]  # keep only relevant columns
    df_filtered.to_csv(outname, index=False)
    print(f"Source data saved to {outname}")


def compute_longitudinal_icc(df_pt, rois_cortical_centiles, outname):
    '''
    Compute ICC for longitudinal data.
    : param df_pt: DataFrame containing data for two timepoints
    : param rois_cortical_centiles: List of ROIs to compute ICC for (e.g., ['centile_CT_bankssts'])
    : param outname: Output filename for the ICC results
    
    : return: DataFrame with ICC values for each ROI
    '''
    
    icc_results = {}

    for region in df_pt[rois_cortical_centiles]:
        if region == 'centile_SA_entorhinal':
            continue
        
        # compute ICC
        icc = pingouin.intraclass_corr(
            data=df_pt, 
            targets='participant', 
            raters='timepoint', 
            ratings=region, 
            nan_policy='omit')
        
        # Extract ICC(3,1): timepoint treated as fixed effect
        icc_row = icc.loc[(icc['Type'] == 'ICC3') & (icc['CI95%'].notnull())].iloc[0]
        icc_value = icc_row['ICC']
        p_value = icc_row['pval']
        ci_lower, ci_upper = icc_row['CI95%']
        df1 = icc_row['df1']
        df2 = icc_row['df2']

        # store into dictionary
        icc_results[region] = (icc_value, df1, df2, ci_lower, ci_upper, p_value)

    # Convert to df
    icc_df = pd.DataFrame.from_dict(icc_results, orient='index', columns=['ICC', 'DOF1', 'DOF2', 'CI95%_lower', 'CI95%_upper', 'pvalue'])
    icc_df['p_fdr'] = multipletests(icc_df['pvalue'], method='fdr_bh')[1]
    
    # format p-values in scientific notation x.xx × 10^-x
    icc_df['ICC'] = icc_df['ICC'].apply(lambda x: f"{x:.3f}")
    icc_df['pvalue'] = icc_df['pvalue'].apply(lambda p: f"{p:.2e}")
    icc_df['p_fdr'] = icc_df['p_fdr'].apply(lambda p: f"{p:.2e}")

    # save to csv
    icc_df.to_csv(outname, index=True)

    return icc_df


def plot_icc_surf(icc_df, outname, brain_measure='SA', limits=(0.5, 1)):
    # get ICC values
    icc_df_val = icc_df['ICC'].to_numpy()
    # insert 0 for entorhinal cortex if SA
    if brain_measure == 'SA':
        icc_df_val = np.insert(icc_df_val, 4, 0)  # insert nan for entorhinal
    # duplicate values for left and right hemisphere
    icc_dupl = np.tile(icc_df_val, (1,2)).squeeze()

    # plot brain map
    plot_brain_map(icc_dupl, outname, cmap='YlOrBr', limits=limits, scale=10, fill=0)