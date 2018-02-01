#!/usr/bin/bash

SCENE_FILENAME="effects_cohens_d_whs26.scene"
WHSIM="26"
IDX=(1 2 3) 
THRESH_NAMES=("small" "medium" "large") 

unset CONTRASTS 
while IFS=$'\r\n' read VALUE; do
    CONTRASTS+=($VALUE)
done < effect_names_whs26.txt

# Change scene file and print
ITER=1
for i in ${IDX[@]}; do 
    
    THRESH=${THRESH_NAMES[$i]} 

    for c in ${CONTRASTS[@]}; do 

        echo "printing $THRESH $c" 
        #wb_command -show-scene effects_cohens_d_whs26.scene $ITER \
        wb_command -show-scene effects_cohens_d_whs26.scene ${THRESH}_${c} \
           ../../images/wb_effects/${c}_effect_cohens_d_whs26_${THRESH}.png \
            10 10 -use-window-size  
    done 
    ITER=$(expr $ITER + 1)
done  

