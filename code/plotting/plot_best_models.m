function [] = plot_best_models( fig, models_to_plot, bestmodel, name, descriptions, R, NS, plotby, filename)

% This function creates area plots of model preference. 
% The plots can have either: 
%   plotby = 'roi': we visualize the percent of subjects that prefer a
%       model for each roi
%   plotby = 'subj': we visualize the percent of rois that prefer a
%       model for each subject

% INPUT: 
%   fig: figure number 
%   models_to_plot: models we want to plot
%   bestmodel: matrix listing which models are best for each subject and roi
%   name: title of the figure window: 
%   descriptions: model descriptoins
%   R: the number of rois
%   NS: the number of subjects
%   plotby: controls whether roi or subject is on the x-axis
%   filename: to save the figure by 

%% Create Area Plot 
figure( fig ); clf;
if strcmp(plotby, 'roi')
    py = zeros( R , length(models_to_plot) );
    for j=1:R
        ps = hist( bestmodel(:,j) , models_to_plot );
        py( j , : ) = ps * 100 / sum( ps );
    end
    % Rearranging ROIs according to the areas best fit by first model
    [ ~ , index ] = sort( py(:,1) , 1 , 'descend' );
    area( py( index , : ) );
    legend( descriptions );
    
    
    xlabel( 'ROI' );
    ylabel( '% Subjects' );
    ylim( [ 0 100 ] );
    xlim( [ 1 R+0.5 ] );
    
elseif strcmp(plotby, 'subj')
    py = zeros( NS , length(models_to_plot) );
    for j=1:NS
        ps = hist( bestmodel(j,:) , models_to_plot );
        py( j , : ) = ps * 100 / sum( ps );
    end
    % Rearranging ROIs according to the areas best fit by first model
    [ ~ , index ] = sort( py(:,1) , 1 , 'descend' );
    area( py( index , : ) );
    legend( descriptions(models_to_plot) );
    xlabel( 'Subject' );
    ylabel( '% ROIs' );
    ylim( [ 0 100 ] );
    xlim( [ 1 NS ] );
end

title( sprintf('Best Model for each ROI and Subject'));
set( gcf , 'Name' , name );
set( gca, 'FontSize', 13);

 print( filename, '-depsc');