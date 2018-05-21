#!/usr/bin/bash

SCENE_FILENAME="cohens_d_whs26_mdl2.scene"
WHSIM="26"
IDX=1
THRESH=0.2 
THRESH_NAMES="small"

unset CONTRASTS 
while IFS=$'\r\n' read VALUE; do
    CONTRASTS+=($VALUE)
done < contrast_names_whs26_mdl2.txt

# Change scene file and print
## change data file thresholds 
source workbench_contrasts_format_mdl2.bash 0.2
THRESH_NAME="small"

# print images from scene file 
ITER=1
for c in ${CONTRASTS[@]}; do 
    wb_command -show-scene cohens_d_whs26_mdl2_small.scene $c ../../images/wb_cohensd/${c}_cohens_d_whs26_mdl2_${THRESH_NAME}.png \
        10 10 -use-window-size  
    ITER=$(expr $ITER + 1)
done 

