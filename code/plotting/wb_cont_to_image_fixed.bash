#!/usr/bin/bash

WHSIM="26"
IDX=( 1 2 3 )
THRESH=(0.2 0.5 0.8) 
THRESH_NAMES=("small" "medium" "large") 

unset CONTRASTS 
while IFS=$'\r\n' read VALUE; do
    CONTRASTS+=($VALUE)
done < contrast_names_whs26.txt

# Change scene file and print
for i in ${IDX[@]}; do 

     echo "$i $THRESH[$i] ${THRESH_NAMES[$i]}"
     ## change data file thresholds 
     source wb_cont_format_fixed.bash ${THRESH[$i]}
     THRESH_NAME=${THRESH_NAMES[$i]}
 
     # print images from scene file 
     ITER=1
     for c in ${CONTRASTS[@]}; do 
         wb_command -show-scene cohens_d_whs26_fixed_${THRESH_NAME}.scene $c ../../images/wb_cohensd/${c}_cohens_d_whs26_fixed_${THRESH_NAME}.png \
             10 10 -use-window-size  
         ITER=$(expr $ITER + 1)
     done 
done  

