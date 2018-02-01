#!/usr/bin/bash 

# RUN from directory ROI2NIFTI/files

# OPTIONS 
SPEC_FILENAME="effects.spec"
WHSIM="26"
FILENAMES=("../../ROI2NIfTI/files/effect_cohensd_small_whs26_whmodel1.dscalar.nii" \
    "../../ROI2NIfTI/files/effect_cohensd_medium_whs26_whmodel1.dscalar.nii" \
    "../../ROI2NIfTI/files/effect_cohensd_large_whs26_whmodel1.dscalar.nii")

# PALETTE 
PALETTE_MODE="MODE_USER_SCALE"
PALETTE_NAME="JET256"
POS_MIN_USER=1.00
POS_MAX_USER=3.00
NEG_MIN_USER=0.00
NEG_MAX_USER=0.00
DISP_NEG=false
DISP_POS=true
DISP_ZERO=0
THRESH_TYPE="THRESHOLD_TYPE_OFF"
THRESH_TEST="THRESHOLD_TEST_SHOW_OUTSIDE"

# Read in Contrast Names 
unset CONTRASTS
while IFS=$'\r\n' read CONTRAST; do
    CONTRASTS+=($CONTRAST)
done < effect_names_whs26.txt

echo "Contrasts:"
for CONTRAST in "${CONTRASTS[@]}"
do
    echo "\t$CONTRAST"
done

# clear file 
rm $SPEC_FILENAME

# add brain files
echo "Adding brain files" 
wb_command -add-to-spec-file $SPEC_FILENAME \
    CORTEX_RIGHT ../../ROI2NIfTI/S1200.R.pial_MSMAll.32k_fs_LR.surf.gii
wb_command -add-to-spec-file $SPEC_FILENAME \
    CORTEX_LEFT ../../ROI2NIfTI/S1200.L.pial_MSMAll.32k_fs_LR.surf.gii

# add all results files 
ITER=1
for FILENAME in ${FILENAMES[@]}; do
    echo "Adding file $FILENAME" 
    wb_command -add-to-spec-file $SPEC_FILENAME CORTEX_LEFT $FILENAME

    echo "\tFormatting" 
    wb_command -cifti-palette $FILENAME $PALETTE_MODE $FILENAME \
        -pos-user $POS_MIN_USER $POS_MAX_USER \
        -neg-user $NEG_MIN_USER $NEG_MAX_USER \
        -palette-name $PALETTE_NAME \
        -disp-pos $DISP_POS \
        -disp-neg $DISP_NEG \

    echo "\tSetting map names" 
    ITER2=1
    for CONTRAST in ${CONTRASTS[@]}
    do  
        wb_command -set-map-name $FILENAME $ITER2 $CONTRAST 
        ITER2=$(expr $ITER2 + 1)
    done
    ITER=$(expr $ITER + 1)
done

# print file information 
# wb_command -file-information $FILENAME

