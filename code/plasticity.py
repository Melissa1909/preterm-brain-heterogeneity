import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from analysis_helpers import get_roi_names


def conduct_pca(data, brain_measure='CT'):
    '''
    Conduct PCA with 3 components on centile scores of cortical ROIs.
    
    data: pd.DataFrame, centile scores for each subject and ROI
    brain_measure: str, 'CT' or 'SA'
    
    returns: np.array, PC1, PC2, PC3 scores for each subject
    '''
    # get ROI names
    rois_cortical = get_roi_names(brain_measure)
    
    # filter df for cortical ROIs
    centile_cols = ['centile_' + roi for roi in rois_cortical]
    data_stripped = data[centile_cols].dropna(axis=1)  # drop SA_entorhinal if SA
    
    # conduct PCA
    print(f'Conducting PCA on {brain_measure} deviation scores...')
    pca = PCA(n_components=3, random_state=123)
    pc_dev_scores = pca.fit_transform(data_stripped.to_numpy())
    
    # explained variance
    explained_variance = pca.explained_variance_ratio_ * 100 
    print(f'Explained variance per principal component: {explained_variance[0]:.2f}%, {explained_variance[1]:.2f}% and {explained_variance[2]:.2f}%')
    
    # pca loadings
    loadings_pc1 = pca.components_[0]
    loadings_pc1 = pd.Series(loadings_pc1, index=data_stripped.columns, name='PC1_loadings')
    
    return pc_dev_scores, loadings_pc1


def plot_pc_loadings(loadings, out_dir, brain_measure_name, ylim=(0, 0.25)):
    '''
    Plot the loadings of PC1 per ROI.

    loadings: df with loadings per ROI
    out_dir: output directory
    brain_measure: str, 'CTh' or 'SA'
    ylim: tuple, y-axis limits
    '''

    fig, ax = plt.subplots(figsize=(8,4))
    sns.barplot(x='ROI', y='PC1_loadings', data=loadings, hue='PC1_loadings', dodge=False, palette='RdYlBu_r', ax=ax)

    ax.get_legend().remove()
    
    # labels
    ax.set_xlabel('')
    xticklabels = [x.get_text().split('_')[2] for x in ax.get_xticklabels()]  # remove centile_CTh from x-axis
    ax.set_xticklabels(xticklabels, rotation=90, fontsize=10)
    ax.set_ylabel('PC loading', fontsize=12)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=10)
    ax.set_ylim(ylim)
    
    ax.set_title(f'{brain_measure_name} PC1 loadings per region', fontsize=14)

    sns.despine()

    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, f'{brain_measure_name}_PC1_loadings.svg'), format='svg', dpi=300)


def read_process_output(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Step 2: Find the start of the table
    start_index = None
    end_index = None

    for i, line in enumerate(lines):
        if "Data for visualizing the conditional effect of the focal predictor:" in line:
            start_index = i + 2  # Table starts two lines below the header
        if "********************" in line and start_index is not None:
            end_index = i
            break

    # Step 3: Extract the table data
    table_lines = lines[start_index:end_index]
    table_data = [line.strip().split() for line in table_lines]

    # Step 4: Convert to a DataFrame
    columns = ["GA", "SES_at_birth", "PC1_CT"]
    df = pd.DataFrame(table_data, columns=columns)

    # Convert columns to numeric types
    df = df.apply(pd.to_numeric, errors='coerce').dropna(how='all')

    return df