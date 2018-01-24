function [] = plot_predicted_vs_actual( tcn, fig, name, models, bestmodel, ...
    models_to_plot, design, var_method, filenm, sub_nums )

% This function plots model predictions against the actual data. It groups
% data points (e.g., a time series from a subject and roi) by which model
% they prefer. We plot predictions for the mean and variance of the average 
% time series at each point. 

%Explanation in paper: 
% We perform predictive checks to qualitatively capture which models 
% accurately describe the mean and variance in the data. We compare three 
% types of time series: 1) model predictions, 2) data generated from the 
% model, and 3) the raw data. We group time series by whether they prefer 
% the \textit{Var+Mean} model or another model. Next, we compute a group 
% mean time series. For each region, we average the time series for all 
% subjects within a group. Then we average each mean regional time series 
% to obtain an overall mean time series. We compute a group mean predicted
% time series using the same procedure but using model predictions instead 
% of the raw data. The procedure changes to when computing the variance. 
% We make the assumption that the set of time series that prefer a given 
% model are samples from a true underlying model, i.e., that we can compute
% variance across regions and subjects rather than only for a single region
% from a single subject. This assumption allows us to estimate the variance 
% at each point in the time series. For each region, we compute the variance 
% over subjects in a group at each time point. Then we compute the mean 
% variance over regions. We plot the actual and predicted mean and variance
% of each of the three time series we compare. The plots serve as a sanity
% check to make sure 1) the model can accurately predict that data and 
% 2) that subjects and regions that prefer the Var+Mean model exhibit 
% non-constant patterns in the variance and 3) that the model can predict 
% these trends.

% INPUT: 
%   tcn: data
%   fig: figure number
%   name: name of figure window 
%   models: model results 
%   bestmodel: matrix of which model is best for each data point 
%   models_to_plot: which models to plot, 
%   design: design matrix with all columns (mean and var) 
%   var_method: method for computing variance
%   filenm: for saving 

%% Plot predicted versus actual 
h = figure( fig ); clf;
h.Position =  [137 47 947 592];
hs = tight_subplot(1,2,[.05 .03],[0.03 .05],[.13 .01]);
set( gcf , 'Name' , name );
lw = 1.5; 
T = size(design, 1); 
models = models(models_to_plot); 
M = length(models); 
NS = size(models{models_to_plot(1)}.allparams, 1); 
R = size(models{models_to_plot(1)}.allparams, 3); 

meandesignnow = design(:, models{1}.meancols);
vardesignnow = design(:, models{1}.varcols);

% plot mean design
for c = 1:size(meandesignnow,2)
    axes(hs(1));
    plot( 1:T , meandesignnow(:,c) * 0.2 + c - 1, '-' , 'LineWidth' , lw ); hold on;
    
end

% plot variance design
for c = 1:size(vardesignnow,2)
    axes(hs(2));
    plot( 1:T , vardesignnow(:,c) * 0.2 + c - 1, '-' , 'LineWidth' , lw ); hold on;
end

% number of conditions to plot
NC = max( size(meandesignnow,2), size(vardesignnow,2));

% get design labels
labels = flip(models{1}.designlabels(1:end));
labels = cellfun( @(x) strrep( x, '_', 'emp '), labels, 'UniformOutput', 0);

for k=1:M
    
    m = models_to_plot(k); 
    allpredsm = models{m}.allpredsm;
    allpredsv = models{m}.allpredsv;
    
    % Find the ROIs and Subject combinations for which model m is best
    wh = find( bestmodel(:) == m );
    % What is the percentage of cases?
    pc = length( wh ) * 100 / ( NS * R );
    titlestr  = sprintf( '%s (%3.2f%%)' , models{ m }.description{1} , pc );
    labels = [labels titlestr];
    
    % Calculate mean and std of data for these cases (separately for each time point)
    meany = nanmean( tcn(:,wh) , 2 );
    % stdy = std( tcn(:,wh) , [] , 2 );
    stdy = var( tcn(:,wh) , [] , 2 );
    
    % Predictions from model for these cases
    ypreds_mean = nanmean( models{ m }.allpredsm(:,wh) , 2 );
    
    %if strcmp( var_method, 'meanpred')
        ypreds_std  = nanmean( models{ m }.allpredsv(:,wh) , 2 );
    %elseif strcmp( var_method, 'sample')
        samples = zeros( T, NS, R);
        for s = 1:NS
            for r = 1:R
                samples(:,s,r) = normrnd( allpredsm(:,s,r), allpredsv(:,s,r));
            end
        end
        samp_ypreds_std = var(samples(:,wh), [], 2);
        samp_ypreds_mean = nanmean( samples(:,wh), 2 );
    %end
    
    % Plot the means
    axes(hs(1));
    plot( 1:T , ypreds_mean - m, 'b-' , 'LineWidth' , lw ); hold on;
    plot( 1:T , samp_ypreds_mean - m, 'b--' , 'LineWidth' , lw ); hold on;
    plot( 1:T , meany - m, 'r-' , 'LineWidth' , lw );
    title( titlestr );
    xlim( [ 0 T ] );
    ylim( [-M-1 NC]);
    
    % Plot the stds
    axes(hs(2));
    plot( 1:T , ypreds_std - m - 1, 'b-' , 'LineWidth' , lw ); hold on;
    plot( 1:T , samp_ypreds_std - m - 1, 'b--' , 'LineWidth' , lw ); hold on;
    plot( 1:T , stdy - m - 1, 'r-' , 'LineWidth' , lw );
    xlim( [ 0 T ] );
    ylim( [-M-1 NC]);
end

labels = flip(labels);

axes(hs(1));
set(gca, 'YTick', -M:(NC-1));
set(gca, 'YTickLabel', labels);
set(gca, 'XTick', []);
set(gca, 'FontSize', 12);
title( 'Mean' );

axes(hs(2));
set(gca, 'YTick', []);
set(gca, 'XTick', []);
set(gca, 'FontSize', 12);
title( 'Variance' );

print(filenm, '-depsc');