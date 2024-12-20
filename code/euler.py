import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



def euler_hemi_combination(euler_file, outdir, dataset_name):
    '''
    Combine Euler indices across hemispheres and save to a new file.
    
    euler_file: file name containing Euler indices for each hemisphere, extracted with mris_euler_number from FreeSurfer
    outdir: os.path, directory to save the combined Euler indices
    dataset_name: str, name of the dataset (e.g., 'BLS-26')
    '''
    # get the sum of Euler characteristics (left hemisphere + right hemisphere)
    subject_sums = {}
    failed_subs = []
    with open(euler_file, "r") as file:
        for line in file:
            parts = line.split()
            subject_id = parts[0]  # extract subject-id
            try:
                value = int(parts[2])  # extract euler number
            except:
                # some subjects have no Euler index as Freesurfer segmentation seems to have failed
                failed_subs.append(subject_id)
                continue
                
            
            # create a dictionary with subject_id: sum of Euler numbers
            if subject_id in subject_sums:
                # Add the value to the existing sum
                subject_sums[subject_id] += value
            else:
                # create a new entry for the subject
                subject_sums[subject_id] = value

    # save this sum across hemis to a new file
    combined_euler_file = os.path.join(outdir, f'QC_output_combined_{dataset_name}.txt')
    with open(combined_euler_file, "w") as file:
        for subject_id, total_sum in subject_sums.items():
            file.write(f"{subject_id} {total_sum}\n")


def get_euler_outliers(combined_euler_file, outdir, threshold=1.5):
    '''
    Determine outliers based on Euler index and threshold*MAD criterion. Will also plot a histogram of Euler indices.
    
    combined_euler_file: file with Euler indices combined across hemispheres
    outdir: directory to save the list of outliers and euler index histogram
    threshold: threshold for exclusion (e.g., 1.5), default is 1.5. Subjects with Euler index greater than threshold*MAD are excluded.
    '''
    # load the euler file with combined hemispheres
    df = pd.read_csv(combined_euler_file, sep=" ", header=None, names=["participant", "euler_index"])

    # calculate median of dataset
    median = df['euler_index'].median()
    print("Median of Euler indices:", median)

    # calculate the threshold for exclusions
    mad = np.median(np.abs(df['euler_index'] - median))
    print("Median absolute deviation (MAD):", mad)
    threshold = threshold * mad

    # filter out subjects with Euler index greater than the threshold
    outlier_df = df[np.abs(df['euler_index'] - median) > threshold]
    outliers = outlier_df['participant'].to_list()
    #  print("Outliers based on Euler index:", outliers)

    # save
    with open(os.path.join(outdir, "euler_outliers.txt"), 'w') as file:
        for subject in outliers:
            file.write(str(subject) + '\n')

    # plot the distribution of Euler indices
    plt.hist(x=df['euler_index'], bins=60)
    plt.axvline(x=median-threshold, color='red', linestyle='--', linewidth=1)
    plt.title('Distribution of Euler Indices')
    plt.xlabel('Euler Index')
    plt.ylabel('Number of Subjects')
    outname=os.path.join(outdir, "euler_histogram.png")
    plt.savefig(outname)
    plt.close()
    print("Euler index histogram saved to:", outname)
    
    return outliers