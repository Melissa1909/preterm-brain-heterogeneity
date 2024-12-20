import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


import sys
sys.path.append('../code/')
from plotting import plot_individual_brain_maps, prepare_data_for_brain_plot


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


# def plot_raw_longitudinal_data_subplots(df, rois, modality, out_dir_fig=os.getcwd()):
#     '''
#     Plot longitudinal data for raw CT or SA in the rois list in individual subplots.
    
#     df: DataFrame containing data for two timepoints, centiles and raw CT
#     rois: List of ROIs to plot (should be 34 in this case)
#     modality: Modality of the data (e.g., CT)
#     out_dir_fig: Output directory to save the figure
#     '''
#     # Define plot layout
#     num_plots = 34
#     num_cols = 7
#     num_rows = (num_plots + num_cols - 1) // num_cols 

#     fig, axes = plt.subplots(num_rows, num_cols, figsize=(25, num_rows*5))
#     axes = axes.flatten()

#     for i, (roi, ax) in enumerate(zip(rois, axes)):
#         # Define centile column based on the ROI
#         centile_column = 'centile_' + roi
#         #centile_column = roi

#         # Filter subjects with extranormal centile scores at any timepoint
#         extranormal_subjects = df.groupby('participant').filter(
#             lambda x: (x[centile_column] < 0.05).any() or (x[centile_column] > 0.95).any()
#         )

#         for sub in extranormal_subjects['participant'].unique():
#             sub_data = extranormal_subjects[extranormal_subjects['participant'] == sub]
            
#             # Iterate through each timepoint for this subject
#             for _, row in sub_data.iterrows():
#                 timepoint = row['timepoint']
#                 ct_value = row[roi]
#                 centile_value = row[centile_column]

#                 # Set color based on centile score
#                 if centile_value < 0.05:
#                     color = 'blue'
#                 elif centile_value > 0.95:
#                     color = 'red'
#                 else:
#                     color = 'k'  # default to black if not extranormal

#                 # Plot raw CT value for this timepoint with appropriate color
#                 ax.scatter(x=timepoint, y=ct_value, color=color, s=40, alpha=.8)

#             # Connect the points with a line
#             ax.plot(sub_data['timepoint'], sub_data[roi], alpha=.8, color='k')

#         # Handle missing values (NaN) and only set y-limits if valid data exists
#         if extranormal_subjects[roi].notna().any():
#             ax.set_ylim(extranormal_subjects[roi].min() * 0.8, extranormal_subjects[roi].max() * 1.2)
#         ax.grid(False)
#         ax.set_xticks([1, 2])
#         ax.tick_params(axis='both', which='major', labelsize=12)
        

#     # Remove unused subplots, if any
#     for j in range(i + 1, len(axes)):
#         fig.delaxes(axes[j])

#     # Save the figure
#     filename = f'{modality}_longitudinal_raw_extranormal.svg'
#     plt.tight_layout()
#     plt.savefig(os.path.join(out_dir_fig, filename))
#     plt.show()



# def split_timepoints(df, roi):  
#     '''
#     Split df into two dataframes based on the timepoint column and merge them on the participant column.
#     df: dataframe containing data for two timepoints
#     roi: the name of the ROI to split and merge
#     '''
#     # split the data by timepoints
#     timepoint1 = df[df['timepoint'] == 1][['participant', roi]]
#     timepoint2 = df[df['timepoint'] == 2][['participant', roi]]

#     # Merge the data on 'participant' to compare the same person's values at both timepoints
#     merged = pd.merge(timepoint1, timepoint2, on='participant', suffixes=('_t1', '_t2'))
#     return merged
    


# def count_deviations_staying_extra(df, rois):
#     '''
#     Count the number of subjects that have infra- or supranormal values at both timepoints.
#     df: dataframe containing data for two timepoints
#     rois: list of ROIs to analyze
#     '''
#     result = []
#     for roi in rois:
#         merged = split_timepoints(df, roi)
#         # Apply conditions for infra and supra significant deviations
#         significant_t1_infra = merged[f'{roi}_t1'] < 0.05
#         significant_t1_supra = merged[f'{roi}_t1'] > 0.95
#         significant_t2_infra = merged[f'{roi}_t2'] < 0.05
#         significant_t2_supra = merged[f'{roi}_t2'] > 0.95
#         # Count how many subjects meet the conditions
#         count_infra_both = (significant_t1_infra & significant_t2_infra).sum()
#         count_supra_both = (significant_t1_supra & significant_t2_supra).sum()
        
#         # Calculate percentages
#         pct_infra_both = 100 if count_infra_both == 0 and significant_t1_infra.sum() == 0 else (count_infra_both / significant_t1_infra.sum()) * 100
#         pct_supra_both = 100 if count_supra_both == 0 and significant_t1_supra.sum() == 0 else (count_supra_both / significant_t1_supra.sum()) * 100
        
#         row = {
#             'ROI': roi,
#             'n_infra_t1': significant_t1_infra.sum(),
#             'n_supra_t1': significant_t1_supra.sum(),
#             'n_infra_t2': significant_t2_infra.sum(),
#             'n_supra_t2': significant_t2_supra.sum(),
#             'n_both_infra': count_infra_both,
#             'n_both_supra': count_supra_both,
#             'pct_infra_both': pct_infra_both,
#             'pct_supra_both': pct_supra_both
#         }
#         result.append(row)
    
#     return pd.DataFrame(result)


# def plot_longitudinal_brain_maps(sub_id, measure, df, outdir_time_1, outdir_time_2):
#     '''
#     Plot individual brain maps showing centile scores for two timepoints.

#     sub_id: subject id, e.g., BEST-BN-001
#     measure: measure to be plotted, e.g., CT
#     df: dataframe containing data for two timepoints
#     outdir_time_1: directory where figures for timepoint 1 should be saved to
#     outdir_time_2: directory where figures for timepoint 2 should be saved to
#     '''
#     # create output directories if they do not exist
#     os.makedirs(outdir_time_1, exist_ok=True)
#     os.makedirs(outdir_time_2, exist_ok=True)

#     # filter for sub_id
#     #df = df[df['participant'] == sub_id]

#     # plot individual brain maps for each timepoint
#     print(f'Plotting individual brain maps for subject {sub_id} timepoint 1...')
#     data_time_1 = df[df['timepoint'] == 1]
#     plot_individual_brain_maps(data_time_1, sub_id, outdir_time_1, scale=10, cmap='RdBu_r', limits=(0, 1), infra_supra=True, measure=measure)
    
#     print(f'Plotting individual brain maps for subject {sub_id} timepoint 2...')
#     data_time_2 = df[df['timepoint'] == 2]
#     plot_individual_brain_maps(data_time_2, sub_id, outdir_time_2, scale=10, cmap='RdBu_r', limits=(0, 1), infra_supra=True, measure=measure)



