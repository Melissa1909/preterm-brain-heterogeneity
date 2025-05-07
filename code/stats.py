import os
import pandas as pd
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