#!/usr/bin/bash 

WHSIM=$1 
DIR=$2

scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_all_subjs.mat \
    ../../Results/${DIR}whs${WHSIM}_all_subjs.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_allllsm.mat \
    ../../Results/${DIR}whs${WHSIM}_allllsm.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_allmodels_2.mat \
    ../../Results/${DIR}whs${WHSIM}_allmodels_2.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_allmodels_4.mat \
    ../../Results/${DIR}whs${WHSIM}_allmodels_4.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_bestmodelCV.mat \
    ../../Results/${DIR}whs${WHSIM}_bestmodelCV.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_allbicm.mat \
    ../../Results/${DIR}whs${WHSIM}_allbicm.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_allmodels_1.mat \
    ../../Results/${DIR}whs${WHSIM}_allmodels_1.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_allmodels_3.mat \
    ../../Results/${DIR}whs${WHSIM}_allmodels_3.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_bestmodelBIC.mat \
    ../../Results/${DIR}whs${WHSIM}_bestmodelBIC.mat
scp ggaut@hpc.oit.uci.edu:/pub/ggaut/VDGLM/Results/${DIR}/whs${WHSIM}_sub_nums.mat \
    ../../Results/${DIR}whs${WHSIM}_sub_nums.mat

