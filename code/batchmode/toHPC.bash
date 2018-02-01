#!/usr/bin/bash

# send main functions to the HPC 
scp analyzedata_batch_v2.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/analyzedata_batch_v2.m
scp null_sample_hypothesis_test.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/null_sample_hypothesis_test.m
scp solve_glm.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/solve_glm.m
scp fit_models_cv.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/fit_models_cv.m
scp fit_models_cv_null_sample.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/fit_models_cv_null_sample.m
scp combine_results.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/combine_results.m
scp combine_results_null_sample_hypothesis_test.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/combine_results_null_sample_hypothesis_test.m
scp find_jobs_to_run.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/find_jobs_to_run.m
scp load_data_and_design.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/load_data_and_design.m
scp set_analysis_options_v2.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/set_analysis_options_v2.m

# optimization functions 
cd ../optimization 
scp loglik_varmean_matrix_var.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/loglik_varmean_matrix_var.m
scp preds_varmean_matrix_var.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/preds_varmean_matrix_var.m
scp loglik_varmean_matrix_logtransform.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/loglik_varmean_matrix_logtransform.m
scp varconstraint2.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/varconstraint2.m

# fmri utility functions
cd ../Utils
scp autocorr_woolrich.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/autocorr_woolrich.m
scp spm_Gpdf.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/spm_Gpdf.m
scp spm_hrf.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/spm_hrf.m
scp set_results_directory.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/set_results_directory.m
cd ../batchmode 

# matlab utility functions 
scp ~/Dropbox/MATLAButils/log4m.m ggaut@hpc.oit.uci.edu:/data/users/ggaut/VDGLM/log4m.m

