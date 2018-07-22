function map = fsl_colormap(varargin)

% This function computes the matlab colors for the FSLview map
% defined in the Human Connectome Project Workbench Viewer

% Scalar values for each color are taken from the Workbench source code
% at https://github.com/Washington-University/workbench/blob/master/src/Files/PaletteFile.cxx:

% palFSLView.addScalarAndColor( 1.0f, '_red_yellow_interp_yellow');
% palFSLView.addScalarAndColor( 0.00001f, '_red_yellow_interp_red');
% palFSLView.addScalarAndColor( 0.0000099f, '_fslview_zero');
% palFSLView.addScalarAndColor(-0.0000099f, '_fslview_zero');
% palFSLView.addScalarAndColor(-0.00001f, '_blue_lightblue_interp_blue');
% palFSLView.addScalarAndColor(-1.0f, '_blue_lightblue_interp_lightblue');


% this->addColor('_red_yellow_interp_red',  255, 0, 0 );
% this->addColor('_red_yellow_interp_yellow',  255, 255, 0 );
% this->addColor('_blue_lightblue_interp_blue',  0, 0, 255 );
% this->addColor('_blue_lightblue_interp_lightblue',  0, 255, 255 );
% this->addColor('_fslview_zero', 0, 0, 0);

if nargin < 1
    n = size(get(gcf, 'Colormap'), 1);
elseif nargin == 1
    n = varargin{1};
else
    error('function takes 0 or 1 input');
end

%% Define color values and interpolate
scalars = [ 0.00001, 1.0, -0.00001, -1.0, 0.0000099, -0.0000099]';

R = [ 255, 255, 0, 0, 0, 0 ]';
G = [ 0,   255, 0, 255, 0, 0]';
B = [ 0, 0, 255, 255, 0, 0]';

[sort_scalars, ii] = sort(scalars);
values = [R, G, B]/255;
values = values(ii,:);

map = interp1(sort_scalars, values, linspace(-1,1,n), 'linear');
