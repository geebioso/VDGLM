#!/usr/bin/bash 
# RUN from directory ROI2NIFTI/files


# SOURCE SOME GLOBAL VARIABLES 
source wb_global_variables.bash

# OPTIONS 
SPEC_FILENAME="contrasts.spec"
WHSIM="26"
FILENAME="${MAIN_FILE_DIRECTORY}/files/cohensd_whs26_whmodel1.dscalar.nii"

# PALETTE 
PALETTE_MODE="MODE_AUTO_SCALE_PERCENTAGE"
PALETTE_NAME="FSL"
POS_MIN=4.00
POS_MAX=96.00
NEG_MIN=2.00
NEG_MAX=98.00
DISP_NEG=true
DISP_POS=true
DISP_ZERO=0
THRESH_TYPE="THRESHOLD_TYPE_NORMAL"
THRESH_TEST="THRESHOLD_TEST_SHOW_OUTSIDE"
THRESH_MIN=$(echo "-$1" | bc)
THRESH_MAX=$(echo "$1" | bc)

# Read in Contrast Names 
unset VALUES
while IFS=$'\r\n' read VALUE; do
    VALUES+=($VALUE)
done < contrast_names_whs26_mdl1.txt

echo "Contrasts:"
for I in "${VALUES[@]}"
do
    echo "\t$I"
done

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
wb_command -cifti-palette $FILENAME $PALETTE_MODE $FILENAME \
    -pos-percent $POS_MIN $POS_MAX \
    -neg-percent $NEG_MIN $NEG_MAX \
    -palette-name $PALETTE_NAME \
    -disp-pos $DISP_POS \
    -disp-neg $DISP_NEG \
    -thresholding $THRESH_TYPE $THRESH_TEST $THRESH_MIN $THRESH_MAX 
     
# Set map names 
echo "Setting map names" 
ITER=1
for VALUE in ${VALUES[@]}
do  
    wb_command -set-map-name $FILENAME $ITER $VALUE 
    ITER=$(expr $ITER + 1)
done

# print file information 
# wb_command -file-information $FILENAME

