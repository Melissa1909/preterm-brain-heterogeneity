import os 
import pandas as pd

import numpy as np
from abagen import fetch_desikan_killiany
from datetime import datetime


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


def dhcp_strip_id(data):
    ''''
    Longitudinal data are saved in one df with subject_id-ses-id as index. This function strips the ses-id part of the index.
    '''
    subject=[]
    session=[]
    for row in data.index:
        # split the index into subject and session
        split = row.split('_')
        subject.append(split[0])
        session.append(split[1])

    # add to df
    data['Subject_ID'] = subject
    data['Session_ID'] = session

    # adapt dtype
    # data['Session_ID'] = data['Session_ID'].astype(int)

    # reset index
    data.reset_index(inplace=True, drop=True)
    return data


def dhcp_format_interview_date(date_str):
    for fmt in ('%d.%m.%y', '%m/%d/%Y'):
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue
    return pd.NaT


def dhcp_assign_scan_label(row):
    if row['scan_age'] >= 37:
        return 'term_equivalent'
    elif row['scan_age'] < 37:
        return 'fetal'
    else:
        return np.nan


def na():
    return slice(None)


def clean_line(line):
        return [part.strip('"') for part in line.strip().split('\t')]


def dhcp_read_meta(file_path):
    '''
    Function to read in the dHCP data from the 18 months follow-up and other meta information. Data is stored in a txt file 
    and therefore needs to be stripped.
    
    file_path: path to the txt file, such as 'stps01.txt'.
    '''
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