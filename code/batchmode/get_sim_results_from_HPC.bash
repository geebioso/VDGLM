#!/usr/bin/bash 

WHSIM=$1 
DIR=$2

scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_all_subjs.mat whs${WHSIM}_all_subjs.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_allllsm.mat whs${WHSIM}_allllsm.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_allmodels_2.mat whs${WHSIM}_allmodels_2.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_allmodels_4.mat whs${WHSIM}_allmodels_4.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_bestmodelCV.mat whs${WHSIM}_bestmodelCV.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_allbicm.mat whs${WHSIM}_allbicm.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_allmodels_1.mat whs${WHSIM}_allmodels_1.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_allmodels_3.mat whs${WHSIM}_allmodels_3.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_bestmodelBIC.mat whs${WHSIM}_bestmodelBIC.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_sub_nums.mat whs${WHSIM}_sub_nums.mat

