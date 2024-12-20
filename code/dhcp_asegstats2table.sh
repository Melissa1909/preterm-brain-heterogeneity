#!/bin/bash

# copy data into structure that FS asegstats2table understands.

MCRIBS_DIR=$1
DATA_DIR=$2

#-------------------------------------------------------------------
# Step 1: Create subject-session IDs for all subjects
#-------------------------------------------------------------------
while IFS=";" read -r SUBJECT_ID SESSION_ID; do
    echo sub-${SUBJECT_ID}-ses-${SESSION_ID} >> ${MCRIBS_DIR}/subjectlist.txt;
done < <(cut -d ";" -f1,2 $DATA_DIR/dHCP_add_info_data_release2.csv | tail -n +2)

#-------------------------------------------------------------------
# Step 2: Bring data into Freesurfer compatible format
#-------------------------------------------------------------------
# Create a freesurfer directory
mkdir -p ${DATA_DIR}/freesurfer_tmp

# Create a folder for each subject and copy data to stats
for sub in $(cat ${MCRIBS_DIR}/subjectlist.txt); do
    OUT_DIR=${DATA_DIR}/freesurfer_tmp/$sub/stats
    mkdir -p ${OUT_DIR}

    # cp data into sub folder and rename files
    cp ${MCRIBS_DIR}/${sub}_aseg.stats ${OUT_DIR}/aseg.stats
done

#-------------------------------------------------------------------
# Step 3: stats2table 
#-------------------------------------------------------------------
export SUBJECTS_DIR=${DATA_DIR}/freesurfer_tmp

asegstats2table --subjectsfile=${MCRIBS_DIR}/subjectlist.txt -t ${DATA_DIR}/dhcp_aseg_mcribs.txt --skip --all-segs
rm ${MCRIBS_DIR}/subjectlist.txt
rm -rf ${DATA_DIR}/freesurfer_tmp