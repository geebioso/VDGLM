#!/usr/bin/bash

SCENE_FILENAME="hcp.scene"
WHSIM="26"
IDX=1
THRESH=0.2 
THRESH_NAMES="small" 
CONTRASTS=("tfMRI_WM_2BK" "tfMRI_WM_0BK" "tfMRI_WM_2BK-0BK")

# Change scene file and print
for i in ${IDX[@]}; do 

     echo "$i $THRESH[$i] ${THRESH_NAMES[$i]}"
     ## change data file thresholds 
     source wb_hcp_format.bash ${THRESH[$i]}
     THRESH_NAME=${THRESH_NAMES[$i]}
 
     # print images from scene file 
     ITER=1
     for c in ${CONTRASTS[@]}; do 
         wb_command -show-scene hcp.scene ${c} ../../images/wb_cohensd/${c}_cohens_d_hcp_${THRESH_NAME}.png \
             10 10 -use-window-size  
         ITER=$(expr $ITER + 1)
     done 
done  

