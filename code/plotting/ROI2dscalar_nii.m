function ROI2dscalar_nii(vec333, fname, modenow)
% ROI333_to_dscalar_nii(vec333, fname)
%  Put the vector of 333 length into a .dscalar.nii file for visualization.
%  Following two files must be in the same folder as this function:
%   Gordon333.32k_fs_LR.dlabel.nii
%   S1200.MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii

% 171211 by Xiangrui Li 

if numel(vec333) ~= 333, error('Input is not length of 333'); end

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
roi = nii_tool('img', [pth '/Gordon333.32k_fs_LR.dlabel.nii']);
roi = squeeze(roi);
roi = roi(1:59412);


nii = nii_tool('load', [pth '/S1200.MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii']);
% zero out all myelin if we are comparing pct
if strcmp(modenow, 'pct')
    nii.img = zeros(size(nii.img)); 
end
for i = 1:333
    nii.img(1,1,1,1,1,roi==i) = vec333(i);
end

nii_tool('save', nii, fname);
