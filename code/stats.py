import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA

from utils import get_roi_names


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



def median_difference_centiles(group1, group2, n_permutations=1000, random_state=123):
    """
    Perform a permutation test to compare the means of two groups.
    
    Parameters:
    - group1: array-like, first group of data
    - group2: array-like, second group of data
    - n_permutations: int, number of permutations to perform
    
    Returns:
    - p_value: float, p-value from the permutation test
    """
    np.random.seed(random_state)
    observed_diff = np.median(group1) - np.median(group2)
    combined = np.concatenate([group1, group2])
    
    perm_diffs = []
    count = 0
    
    for _ in range(n_permutations):
        np.random.shuffle(combined)
        perm_group1 = combined[:len(group1)]
        perm_group2 = combined[len(group1):]
        perm_diff = np.median(perm_group1) - np.median(perm_group2)
        perm_diffs.append(perm_diff)
        if abs(perm_diff) >= abs(observed_diff):
            count += 1
            
    p_value = count / n_permutations
    
    return p_value, perm_diffs, observed_diff


def plot_permuation_test(observed_diff, perm_diffs, roi, out_dir):
    """
    Plot the distribution of permutation test differences.
    
    Parameters:
    - diffs: array-like, differences from the permutation test
    """
    plt.figure(figsize=(5, 3))
    sns.histplot(perm_diffs, bins=30, kde=True)
    plt.axvline(observed_diff, color='red', linestyle='--', label='Observed median difference')
    plt.title(roi)
    plt.xlabel('Difference in medians (FT - PT)')
    plt.ylabel('Frequency')
    #plt.savefig(os.path.join(out_dir, f'{roi}_median_difference_perm_test.svg'), bbox_inches='tight', dpi=300, format='svg')