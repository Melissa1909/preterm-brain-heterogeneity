import os
from os.path import join
import numpy as np
import pandas as pd

import nibabel as nib
from nibabel.gifti import GiftiImage, GiftiDataArray

from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

from enigmatoolbox.permutation_testing import spin_test
from neuromaps.stats import compare_images
from neuromaps.nulls import burt2020



def get_mean_expression_cell_type(cell_type, expression, cell_type_genes, out_dir):
    '''
    Get mean expression for all genes characteristic for cell_type.
    cell_type: str, cell type (e.g., Astro)
    expression: pd.DataFrame, AHBA gene expression
    cell_type_genes: list, list of genes characteristic for cell_type
    out_dir: os.path, where output should be stored
    '''
    # extract genes for current cell type
    cell_genes = cell_type_genes[cell_type]
    print(f"{len(cell_genes)} genes found for {cell_type}")
    
    # generate empty array
    expressionAll = np.zeros([len(cell_genes), 68])
    used_genes = []
    
    for g, gene in enumerate(cell_genes):
        if gene in expression.index:
            gene_expr = expression.loc[gene].values
            used_genes.append(gene)
        else:
            #print(f'Skipping gene {gene}')
            gene_expr = np.nan

        # store expression values for all genes of relevant cell type
        expressionAll[g,:] = gene_expr
    
    # empty nan rows
    expressionAll = expressionAll[~np.isnan(expressionAll).all(axis=1)]
    expressionAll = pd.DataFrame(expressionAll)
    expressionAll.index = used_genes
    print(f'{len(used_genes)} genes used for {cell_type}')
    
    # save list of used genes
    with open(join(out_dir, f'used_genes_{cell_type}.txt'), 'w') as f:
        for gene in used_genes:
            f.write(f"{gene}\n")
            
    # get mean expression for gene set for each DK-ROI
    expressionAllMean = expressionAll.mean(axis=0)
    
    # save mean expression for cell type
    expressionAllMean_df = pd.DataFrame(expressionAllMean, columns=[cell_type])
    expressionAllMean_df['rois'] = expression.columns
    expressionAllMean_df.to_csv(join(out_dir, f'expression_{cell_type}.csv'), sep=',')
    print(f"Mean expression for {cell_type} saved to {out_dir}/expression_{cell_type}.csv")
    
    return expressionAllMean


def generate_parcellation_file():
    '''
    Generate parcellation file for DK-68 atlas that is required for burt2020.
    
    Attention: requires FreeSurfer to be installed and the environment variable FREESURFER_HOME to be set.
    '''
    
    freesurfer_home = os.environ.get('FREESURFER_HOME')
    
    for hemi in ['lh', 'rh']:
        labels, _, _ = nib.freesurfer.read_annot(
            os.path.join(freesurfer_home, 'subjects', 'fsaverage5', f'label/{hemi}.aparc.annot'), 
            orig_ids=False)
        
        gii_array = GiftiDataArray(data=labels)   
        gii_img = GiftiImage(darrays=[gii_array])
    
        os.makedirs('outputs/tmp', exist_ok=True)
        nib.save(gii_img, f'outputs/tmp/{hemi}.aparc.label.gii')



def cell_correlation(mean_expression, cortical_data, n_rot=1000, p_method='spin'):
    ''''
    Correlate mean gene expression of a specific cell type (mean_expression) with deviation scores for each subject.
    
    mean_expression: np.array, mean expression for cell type from get_mean_expression_cell_type()
    cortical_data: pd.DataFrame (nsub x 68), deviation scores per subject
    n_rot: int, how many times spin_test should be performed
    p_method: which method to use for p-value calculation, 'spin', 'burt', or 'none'
    '''
    # create empty array for storing output
    corr_coeffs = np.zeros([cortical_data.shape[0], 3])
    
    for s, sub in enumerate(cortical_data.index):
        sub_df = cortical_data.loc[sub]
        
        r, p = spearmanr(mean_expression, sub_df)
        
        # calculate p-value with method that corrects for spatial autocorrelation
        if p_method == 'spin':
            pspin = spin_test(mean_expression, sub_df, n_rot=n_rot, type='spearman')
            
        elif p_method == 'burt':
            generate_parcellation_file()
            parcellation = ('outputs/tmp/lh.aparc.label.gii', 'outputs/tmp/rh.aparc.label.gii')
            rotated = burt2020(data=mean_expression, atlas='fsaverage', density='10k',
                                parcellation=parcellation, n_perm=n_rot, seed=1234)
            r, pspin = compare_images(mean_expression, sub_df, nulls=rotated, metric='spearmanr')
            
        # elif p_method == 'burt':
        #     generate_parcellation_file()
        #     parcellation = ('outputs/tmp/lh.aparc.label.gii', 'outputs/tmp/rh.aparc.label.gii')
        #     rotated = burt2020(data=mean_expression, atlas='fsaverage', density='10k',
        #                         parcellation=parcellation, n_perm=n_rot, seed=1234)
        #     r, pspin = compare_images(mean_expression, sub_df, nulls=rotated, metric='spearmanr')
            
        elif p_method == 'none':
            pspin = np.nan
        else:
            raise ValueError("p_method must be 'spin', 'burt', or 'none'")
        
        
        # store output
        corr_coeffs[s,0] = r
        corr_coeffs[s,1] = p
        corr_coeffs[s,2] = pspin
        
        # to df
        corr_coeffs_df = pd.DataFrame(corr_coeffs, columns=['Spearman_r', 'p', f'p{p_method}'], index=cortical_data.index)
        
        # FDR correction
        corr_coeffs_df['p_fdr'] = multipletests(corr_coeffs_df['p'], method='fdr_bh')[1]

    return corr_coeffs_df



