#!/usr/bin/bash

SCENE_FILENAME="tstats_whs26.scene"
WHSIM="26"
IDX=( 1 2 3 )
THRESH_NAMES=("small" "medium" "large") 

unset CONTRASTS 
while IFS=$'\r\n' read VALUE; do
    CONTRASTS+=($VALUE)
done < contrast_names_whs26.txt


## change data file thresholds 

# print images from scene file 
ITER=1
for c in ${CONTRASTS[@]}; do 
    wb_command -show-scene tstats_whs26.scene $c ../../images/wb_tstats/${c}_tstats_whs26.png \
        10 10 -use-window-size  
    ITER=$(expr $ITER + 1)
done 

