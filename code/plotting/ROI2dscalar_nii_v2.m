function ROI2dscalar_nii(rst333, fname, mapName,  modenow, dotest)
% ROI2dscalar_nii(rst333, fname, mapName) % convert one map
% ROI2dscalar_nii(rst333, fname, mapNamesCell) % convert multiple maps
%  Put the result of 333 ROIs (333 by nMap) into a .dscalar.nii file for visualization.
%  fname is the dscalar.nii file to save the result.
%  Optional mapName is used to label the map. If left out, the variable name
%  will used as the map name.
%  Following two files must be in the same folder as this function:
%   Gordon333.32k_fs_LR.dlabel.nii
%   S1200.MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii

% 171211 by Xiangrui Li 
% 180129 Work for multiple maps 

if and( size(rst333, 1) ~= 333, dotest = 0), error('first dim must be 333'); end
nMap = size(rst333, 2);

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

filepath = fullfile( getenv( 'HOME'), 'Dropbox', 'FMRI', 'Projects', 'varianceGLM', 'ROI2NIfTI'); 
roi = nii_tool('img', [filepath '/Gordon333.32k_fs_LR.dlabel.nii']);
roi = squeeze(roi);
roi = roi(1:59412);

nii = nii_tool('load', [filepath '/S1200.MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii']);
nii.img = repmat(nii.img, [1 1 1 1 size(rst333,2) 1]);
if strcmp(modenow, 'pct')
        nii.img= zeros(size(nii.img)); 
end
for i = 1:333
    ind = roi==i;
    img = repmat(rst333(i,:), [sum(ind) 1]);
    nii.img(1,1,1,1,:,ind) = permute(img, [3:6 2 1]);
end



if nargin<3 || isempty(mapName), mapName = inputname(1); end
if isempty(mapName), mapName = 'myMap'; end

% mapName = matlab.lang.makeValidName(mapName);

if ~iscell(mapName), mapName = {mapName}; end
if nMap ~= numel(mapName)
    error('Number of maps in data and mapName don''t match');
end 

str = '';
for i = 1:numel(mapName)
    str = sprintf('%s<NamedMap>\n\t<MapName>%s</MapName>\n</NamedMap>\n', str, mapName{i});
end
[i1, i2] = regexp(nii.ext.edata_decoded, '<NamedMap>.*?</NamedMap>', 'start', 'end');
nii.ext.edata = [nii.ext.edata_decoded(1:i1-1) str nii.ext.edata_decoded(i2+1:end)];

nii.hdr.cal_min = -1; % threshold: negative value evokes two-sided LUT
nii.hdr.cal_max = 3;  % clip value
nii_tool('save', nii, fname);
