#!/usr/bin/bash 

# RUN from directory ROI2NIFTI/files

source wb_global_variables.bash

# OPTIONS 
SPEC_FILENAME="contrasts_fixed_mdl2.spec"
WHSIM="26"
FILENAME="${MAIN_FILE_DIRECTORY}/files/cohensd_whs26_whmodel2.dscalar.nii"

# PALETTE 
PALETTE_MODE="MODE_USER_SCALE"
PALETTE_NAME="FSL"
POS_MIN_USER=0.00
NEG_MAX_USER=0.00
DISP_NEG=true
DISP_POS=true
DISP_ZERO=0
THRESH_TYPE="THRESHOLD_TYPE_NORMAL"
THRESH_TEST="THRESHOLD_TEST_SHOW_OUTSIDE"
THRESH_MIN=$((-$1))
THRESH_MAX=$(($1))

# Read in Contrast Names 
unset CONTRASTS
while IFS=$'\r\n' read CONTRAST; do
    CONTRASTS+=($CONTRAST)
done < contrast_names_whs26_mdl2.txt

echo "Contrasts:"
for CONTRAST in "${CONTRASTS[@]}"
do
    echo "\t$CONTRAST"
done

# Read in Mean Bounds 
unset MEAN_BOUNDS
while IFS=$'\r\n' read VALUE; do
    MEAN_BOUNDS+=($VALUE)
done < bounds_whs26_mdl2_mean.txt

# Read in Var Bounds 
unset VAR_BOUNDS
while IFS=$'\r\n' read VALUE; do
    VAR_BOUNDS+=($VALUE)
done < bounds_whs26_mdl2_var.txt

# clear file 
rm $SPEC_FILENAME

# add brain files
echo "Creating .spec file" 
echo "Adding brain files" 
wb_command -add-to-spec-file $SPEC_FILENAME \
    CORTEX_RIGHT ${MAIN_FILE_DIRECTORY}/S1200.R.pial_MSMAll.32k_fs_LR.surf.gii
wb_command -add-to-spec-file $SPEC_FILENAME \
    CORTEX_LEFT ${MAIN_FILE_DIRECTORY}/S1200.L.pial_MSMAll.32k_fs_LR.surf.gii
wb_command -add-to-spec-file $SPEC_FILENAME \
    CORTEX ${MAIN_FILE_DIRECTORY}/Gordon333.32k_fs_LR.dlabel.nii

# add all results files 
echo "Adding Cohen's d file" 
wb_command -add-to-spec-file $SPEC_FILENAME CORTEX_LEFT $FILENAME

echo "Formatting Cohen's d file" 

# Set map names 
echo "Setting map names" 

ITER=1

for c in ${CONTRASTS[@]}; do

    # Set map name 
    wb_command -set-map-name $FILENAME $ITER $c
    ITER=$(expr $ITER + 1)

    # get bounds for mean or var contrast
    if [[ $c = *"mean"* ]]; then
      POS_MAX_USER=${MEAN_BOUNDS[2]}
      NEG_MIN_USER=${MEAN_BOUNDS[1]}
    else
      POS_MAX_USER=${VAR_BOUNDS[2]}
      NEG_MIN_USER=${VAR_BOUNDS[1]}
    fi

    # format map 
    wb_command -cifti-palette $FILENAME $PALETTE_MODE $FILENAME \
        -column $c \
        -pos-user $POS_MIN_USER $POS_MAX_USER \
        -neg-user $NEG_MIN_USER $NEG_MAX_USER \
        -palette-name $PALETTE_NAME \
        -disp-pos $DISP_POS \
        -disp-neg $DISP_NEG \
        -thresholding $THRESH_TYPE $THRESH_TEST $THRESH_MIN $THRESH_MAX 
done
     
# print file information 
# wb_command -file-information $FILENAME

