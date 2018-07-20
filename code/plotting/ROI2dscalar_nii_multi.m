function ROI2dscalar_nii_multi(cohensd, contrast_names, fname, modenow)
% ROI333_to_dscalar_nii(vec333, fname)
%  Put the vector of 333 length into a .dscalar.nii file for visualization.
%  Following two files must be in the same folder as this function:
%   Gordon333.32k_fs_LR.dlabel.nii
%   S1200.MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii

% 171211 by Xiangrui Li 

% INPUT: 
%   array whsims: numbers of simulations we want to run 
%   cell cohends: contains cohendsd for eaach contrast
%   cell contrast_names: contains the name of each contrast 
%   str fname: filename 
%   str modenow: 

files_pth = fullfile( getenv('HOME'), 'Dropbox', 'FMRI', 'Projects', 'varianceGLM', 'ROI2NIfTI'); 
addpath(fullfile(getenv('HOME'), 'Dropbox', 'FMRI', 'gifti-1.6')); 
addpath(fullfile(getenv('HOME'), 'Dropbox', 'FMRI', 'cifti-matlab')); 

wbcommand =  '/Applications/workbench/bin_macosx64/wb_command'; 

NC = length(cohensd); 

if nargin<2 || isempty(fname)
    [fname, pth] = uiputfile('*.dscalar.nii', 'Input file name to save result');
    if isnumeric(fname), return; end
    fname = fullfile(pth, fname);
end
if isempty(regexp(fname, '.dscalar.nii$', 'once'))
    if ~isempty(regexp(fname, '.nii$', 'once')), fname = fname(1:end-4); end
    if ~isempty(regexp(fname, '.dscalar$', 'once')), fname = fname(1:end-8); end
    fname = [fname '.dscalar.nii'];
end

pth = fileparts(mfilename('fullpath'));

% 
roi = nii_tool('img', [files_pth '/Gordon333.32k_fs_LR.dlabel.nii']);
roi = squeeze(roi);
roi = roi(1:59412);

% nii2 = nii_tool('load', [files_pth '/S1200.MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii']);
% 
cii = ciftiopen([files_pth '/HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s4_MSMSulc.dscalar.nii'],...
     wbcommand); 
% 
% [gifti_obj,xml_data] = cifti_open([pth '/HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s4_MSMSulc.dscalar.nii'],...
%     wbcommand); 

% cii = ft_read_cifti([pth '/HCP_S1200_997_tfMRI_ALLTASKS_level2_cohensd_hp200_s4_MSMSulc.dscalar.nii'])
if strcmp(modenow, 'pct')
        cii.cdata= zeros(size(cii.cdata)); 
end
    
for c = 1:NC 
% zero out all myelin if we are comparing pct
    
    for i = 1:333
        cii.cdata(roi==i, c) = cohensd{c}(i);
    end

end

cii.cdata(:, c+1:end) = []; 

ciftisavereset(cii,fname,wbcommand)
% nii.hdr.dim(6) = size(nii.img, 5); 
% 
% 
% nii_tool('save', nii, fname);
