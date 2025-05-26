# Helpers for preprocessing

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pingouin import ancova
from statsmodels.stats.multitest import multipletests


# def load_freesurfer_aparc(file):
#     '''
#     Load aparc file created by FreeSurfer and combined with aparcstats2table.
    
#     file: directory to the aparc file
#     '''
#     df = pd.read_csv(file, sep='\t')

#     # rename ?h.aparc.modality column to participant
#     df.rename(columns=lambda col: 'participant' if 'aparc' in col else col, inplace=True)
#     df.drop(columns=['BrainSegVolNotVent', 'eTIV'], inplace=True)
#     return df


# def load_freesurfer_aseg(file):
#     '''
#     Load aseg file created by FreeSurfer and combined with asegstats2table.
    
#     file: directory to the aseg file
#     '''
#     df = pd.read_csv(file, sep='\t')

#     # rename Measure:volume column to participant
#     df.rename(columns=lambda col: 'participant' if 'Measure' in col else col, inplace=True)
#     return df


def calculate_scanner_difference(data, rois, scanner_var, covars):
    '''
    Compute the F-values for each ROI to test for scanner differences with covariates in "covars". Applies FDR-correction as well.
    
    data: pd.DataFrame, df containing modality values for all ROIs as well as scanner, age, and binarized sex information
    rois: list of ROI names
    scanner_var: str, between-subject factor (e.g., 'Scanner-ID')
    covars: list of covariates to control for
    '''

    # get F-values and corresponding p-value for each ROI
    f_data_all = np.empty(shape=[0,0])
    p_values_all = np.empty(shape=[0,0])

    for r in rois:
        mdl = ancova(data, dv=r, between=scanner_var, covar=covars)
        f = mdl.to_numpy()[0,3]
        p = mdl.to_numpy()[0,4]
        f_data_all = np.append(f_data_all, f)
        p_values_all = np.append(p_values_all, p)

    # fdr correction
    p_values_fdr = multipletests(p_values_all, alpha=0.05, method='fdr_bh')[1]

    return pd.DataFrame(data = {'ROI': rois, 'F-values': f_data_all, 'p-values': p_values_all, 'p-values_fdr': p_values_fdr})

# WAS NEEDED FOR RELEASE 2 DATA
# def dhcp_aparcstats2table(csv_files, brain_measure):
#     '''
#     Rearange individual freesurfer output files into one table. 
    
#     filename: name of txt file, such as 'stps01.txt'.
#     '''
#     aparc = pd.DataFrame()
    
#     for file in csv_files:
#         # Extract subject name from the file path
#         dirname, subject_name = os.path.split(file)
#         subject_name=subject_name[:-15]
        
#         # Read the stats file
#         df = pd.read_csv(file, sep=',')
        
#         if brain_measure == 'CT':
#             # Set the subject name as the index for the extracted column
#             mri_data_sub = df[["ThickAvg"]].rename(columns = {"ThickAvg": subject_name})



#         elif brain_measure == 'SA':
#             mri_data_sub = df[["SurfArea"]].rename(columns = {"SurfArea": subject_name})
            
#         else:
#             raise ValueError('brain_measure must be either CT or SA')
        
        
#         # Append the extracted column to the summary DataFrame
#         mri_data_sub.index = df['StructName']
#         aparc = pd.concat([aparc, mri_data_sub], axis=1).copy()

#     aparc = aparc.transpose()
#     if 'Medial_Wall' in aparc.columns:
#         aparc.drop('Medial_Wall',axis=1, inplace=True)
        
#     return aparc
    
# WAS NEEDED FOR RELEASE 2 DATA
# def strip_dhcp_id(data):
#     ''''
#     Longitudinal data are saved in one df with subject_id-ses-id as index. This function strips the ses-id part of the index.
#     '''
#     subject=[]
#     session=[]
#     for row in data.index:
#         # split the index into subject and session
#         split = row.split('-')
#         subject.append(split[1])
#         session.append(split[3])

#     # add to df
#     data['Subject_ID'] = subject
#     data['Session_ID'] = session

#     # adapt dtype
#     data['Session_ID'] = data['Session_ID'].astype(int)

#     # reset index
#     data.reset_index(inplace=True, drop=True)
#     return data


def clean_line(line):
        return [part.strip('"') for part in line.strip().split('\t')]


def read_dhcp_18months(filename, meta_dir):
    '''
    Function to read in the dHCP data from the 18 months follow-up. Data is stored in a txt file 
    and therefore needs to be stripped.
    
    filename: name of txt file, such as 'stps01.txt'.
    '''
    file_path=os.path.join(meta_dir, 'dHCP_restrictedLongData','1222825',filename)

    with open(file_path, 'r') as f:
        cleaned_data = [clean_line(line) for line in f]
        
    df = pd.DataFrame(cleaned_data, columns=cleaned_data[0])
    df = df.iloc[2:,:]
    return df


def reorder_vars(first_vars, df, idps):
    '''
    Function from Leon D Lotter to reorder the columns of a dataframe.
    
    first_vars: list of variables that should be at the beginning of the dataframe
    df: dataframe to reorder
    idps: list of independent variables
    '''
    idp_vars = []
    for i in idps:
        idp_vars += [c for c in df.columns if i in c]
    return df[first_vars + [c for c in df.columns if c not in first_vars+idp_vars] + idp_vars].copy()

# def dhcp_aparcstats2table(csv_files, brain_measure):
#     '''
#     Rearange individual freesurfer output files into one table. 
    
#     filename: name of txt file, such as 'stps01.txt'.
#     '''
#     aparc = pd.DataFrame()
    
#     for file in csv_files:
#         # Extract subject name from the file path
#         dirname, subject_name = os.path.split(file)
#         subject_name=subject_name[:-15]
        
#         # Read the stats file
#         df = pd.read_csv(file, sep=',')
        
#         if brain_measure == 'CT':
#             # Set the subject name as the index for the extracted column
#             mri_data_sub = df[["ThickAvg"]].rename(columns = {"ThickAvg": subject_name})

#         elif brain_measure == 'SA':
#             mri_data_sub = df[["SurfArea"]].rename(columns = {"SurfArea": subject_name})
            
#         else:
#             raise ValueError('brain_measure must be either CT or SA')
        
        
#         # Append the extracted column to the summary DataFrame
#         mri_data_sub.index = df['StructName']
#         aparc = pd.concat([aparc, mri_data_sub], axis=1).copy()

#     aparc = aparc.transpose()
#     if 'Medial_Wall' in aparc.columns:
#         aparc.drop('Medial_Wall',axis=1, inplace=True)
        
#     return aparc
    






def na():
    return slice(None)
