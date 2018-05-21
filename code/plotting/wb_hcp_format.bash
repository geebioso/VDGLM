#!/usr/bin/bash

FILENAME="/Users/Garren/Dropbox/FMRI/Projects/varianceGLM/ROI2NIfTI/HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s4_MSMSulc.dscalar.nii"
CORT_R_FILENAME="/Users/Garren/Dropbox/FMRI/Projects/varianceGLM/ROI2NIfTI/S1200.R.pial_MSMAll.32k_fs_LR.surf.gii"
CORT_L_FILENAME="/Users/Garren/Dropbox/FMRI/Projects/varianceGLM/ROI2NIfTI/S1200.L.pial_MSMAll.32k_fs_LR.surf.gii"
SPEC_FILENAME="hcp.spec"

rm $SPEC_FILENAME 

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
THRESH_MIN=-0.2
THRESH_MAX=0.2

wb_command -add-to-spec-file $SPEC_FILENAME \
    CORTEX_RIGHT $CORT_R_FILENAME
wb_command -add-to-spec-file $SPEC_FILENAME \
    CORTEX_LEFT $CORT_L_FILENAME

wb_command -cifti-palette $FILENAME $PALETTE_MODE $FILENAME \
    -pos-percent $POS_MIN $POS_MAX \
    -neg-percent $NEG_MIN $NEG_MAX \
    -palette-name $PALETTE_NAME \
    -disp-pos $DISP_POS \
    -disp-neg $DISP_NEG \
    -thresholding $THRESH_TYPE $THRESH_TEST $THRESH_MIN $THRESH_MAX

wb_command -add-to-spec-file $SPEC_FILENAME CORTEX_LEFT $FILENAME
