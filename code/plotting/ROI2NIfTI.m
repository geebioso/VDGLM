function [h] = ROI2NIfTI(rst, fname, viewer, options)
% ROI2NIfTI(rst, fname)
% 

% GG added 
% viewer: which viewer to use (either 'nii_viewer' or 'view_nii' )


% Convert ROI analysis result into NIfTI for visualization.
% The rst is in dimension of nROIs by nVol, where nVol is often 1.
% 
% fname is the result NIfTI name to save for visualization. If not
% provided, 'myTempRst.nii.gz' will be saved into current directory.
% 
% Example:
%  rst = rand(269, 1); % fake result for 269 ROIs
%  ROI2NIfTI(rst, './myFakeRst.nii.gz') % save to NIfTI
%  % overlay result NIfTI onto standard brain
%  nii_viewer('standard_brain.nii.gz', './myFakeRst.nii.gz')
% 
% Required: NIfTI tool at http://www.mathworks.com/matlabcentral/fileexchange/42997

% 170830 Xiangrui Li

fname = [fname '.gz.nii']; 

addpath('/Users/Garren/Dropbox/FMRI/BrainVisualization/NIfTI/NIfTI');
addpath('/Users/Garren/Dropbox/FMRI/Code/varianceGLM/ROI2NIfTI/dicm2nii'); 

pth = fileparts(mfilename('fullpath')); % path of this code

rst = squeeze(rst);
if size(rst,1)==1, rst = rst'; end
nROI = size(rst, 1);
if nROI == 269 % special case
    load([pth '/ROI_ind_269.mat']); % variable 'roi'
else
    roi = nii_tool('img', [pth '/ROI_299_3mm.nii.gz']); % 299 ROIs
end
ind = unique(roi); ind(1) = [];

nii = nii_tool('load', [pth '/standard_brain_3mm.nii.gz']);
% nii.hdr.intent_code = 3; % 3:T-test; 2: correlation
% nii.hdr.intent_p1 = DoF;
% nii.hdr.descrip = 'any comment';

nii.img(:) = 0; % zero the img
img = nii.img;
nVol = numel(rst) / nROI; % number of volumes 
nii.img = repmat(nii.img, [1 1 1 nVol]);
for j = 1:nVol
    img(:) = 0;
    for i = 1:numel(ind)
        a = roi == ind(i);
        img(a) = rst(i,j);
    end
    nii.img(:,:,:,j) = img;
end

if nargin<2 || isempty(fname), fname = 'myTempRst.nii.gz'; end
nii_tool('save', nii, fname);

switch viewer 
    case 'nii_viewer'
        nii_viewer([pth '/standard_brain_3mm.nii.gz'], fname); 
        h = NaN;
    case 'view_nii'
        nii = load_nii(fname); 
        h = view_nii( nii , options); 
end



