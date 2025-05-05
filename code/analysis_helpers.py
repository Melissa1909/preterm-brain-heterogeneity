# General helpers

import os
import pandas as pd
import math
#import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests
from abagen import fetch_desikan_killiany


def get_roi_names(brain_measure, bilateral=True, global_vars=False):
    '''
    Get the names of the cortical Desikan-Killiany ROIs adding the correct modality and hemisphere prefix if needed.

    brain_measure: brain_measure of the data (i.e., CT, SA)
    bilateral: if True, the hemisphere prefixes will be added to the ROI names
    global_vars: if True, global measures will be added to the output list (i.e., eTIV, WMV, GMV, and sGMV)
    '''

    atlas = fetch_desikan_killiany(surface=True)
    atlas = pd.read_csv(atlas['info'])
    atlas = atlas[(atlas['structure'] == 'cortex') & (atlas['hemisphere'] == 'L')]
    labels = atlas['label'].tolist()

    if bilateral == False:
        rois = ['lh_'+ brain_measure + '_' + label for label in labels] + ['rh_'+ brain_measure + '_' + label for label in labels]
    else:
        rois = [brain_measure + '_' + label for label in labels]

    if global_vars == True:
        rois += ['WMV', 'GMV', 'sGMV']

    return rois


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


