function [] = plot_predicted_vs_actual_varmean_only( tcn, fig, name, models, bestmodel, design, filenm )

% This function is similar to plot_predicted_vs_actual.m except that it
% only plots the Var+Mean model and groups data points by whether a model
% with variance terms is preferred (e.g., Var+Mean or Var) or not. 

h = figure( fig ); clf;
h.Position =  [137 47 947 592];
hs = tight_subplot(1,2,[.05 .03],[0.03 .05],[.13 .01]);
set( gcf , 'Name' , name );
lw = 1.5; 
T = size(design, 1); 
M = length(models); 
NS = size(models{1}.allparams, 1); 
R = size(models{1}.allparams, 3); 

% get the design matrix, this is the same matrix for all
% OSU subjects. For the HCP, there are slight differences
% between subjects, but for simplicity, we plot the design from
% the first subject
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

allpredsm = models{1}.allpredsm;
allpredsv = models{1}.allpredsv;
allpredsv( allpredsv < 0 ) = 0.01;

% Find the ROIs and Subject combinations for which model m is best
% wh = find( bestmodelCV(:) == 1 );
wh = find(ismember( bestmodel(:), [1,3]));
% What is the percentage of cases?
pc = length( wh ) * 100 / ( NS * R );
titlestr  = sprintf( '%s (%3.2f%%)' , models{ 1 }.description{1}, pc );

% Calculate mean and std of data for these cases (separately for each time point)
meany = zeros(T,R);
vary  = zeros(T,R);
samp_ypreds_mean = zeros(T,R);
samp_ypreds_var  = zeros(T,R);
ypreds_mean = zeros(T,R);
ypreds_var  = zeros(T,R);

samples = zeros(T, NS, R);
for r = 1:R
    for s = 1:NS
        samples(:,s,r) = normrnd( allpredsm(:,s,r), sqrt(allpredsv(:,s,r)));
    end
end
for r = 1:R
    wh_subjs = ismember(bestmodel(:,r), [1, 3]);
    meany(:,r)    = nanmean( tcn(:,wh_subjs,r), 2);
    vary(:,r)     = var(  tcn(:,wh_subjs,r), [], 2);
    samp_ypreds_mean(:,r) = nanmean( samples(:, wh_subjs, r), 2);
    samp_ypreds_var(:,r)  = var( samples(:, wh_subjs, r), [], 2);
    
    ypreds_mean(:,r) = nanmean(  models{ 1 }.allpredsm(:, wh_subjs, r), 2);
    ypreds_var(:,r)  = nanmean(  models{ 1 }.allpredsv(:, wh_subjs, r), 2);
end

meany = nanmean( meany, 2 );
vary  = nanmean( vary, 2 );
samp_ypreds_mean = nanmean( samp_ypreds_mean, 2 );
samp_ypreds_var = nanmean( samp_ypreds_var, 2 );
ypreds_mean = nanmean( ypreds_mean, 2 );
ypreds_var = nanmean( ypreds_var, 2 );

labels = [labels 'Var+Mean/Var'];

% Plot the means
axes(hs(1));
plot( 1:T , ypreds_mean - 1, 'b-' , 'LineWidth' , lw ); hold on;
plot( 1:T , samp_ypreds_mean - 1, 'b--' , 'LineWidth' , lw ); hold on;
plot( 1:T , meany - 1, 'r-' , 'LineWidth' , lw );
title( titlestr );
xlim( [ 0 T ] );
ylim( [-M+1 NC]);

% Plot the var
axes(hs(2));
plot( 1:T , ypreds_var - 1 - 1, 'b-' , 'LineWidth' , lw ); hold on;
plot( 1:T , samp_ypreds_var - 1 - 1, 'b--' , 'LineWidth' , lw ); hold on;
plot( 1:T , vary - 1 - 1, 'r-' , 'LineWidth' , lw );
xlim( [ 0 T ] );
ylim( [-M+1 NC]);

% % Models that prefer other models than the V+M model
% Calculate mean and std of data for these cases (separately for each time point)
meany = zeros(T,R);
vary  = zeros(T,R);
samp_ypreds_mean = zeros(T,R);
samp_ypreds_var  = zeros(T,R);
ypreds_mean = zeros(T,R);
ypreds_var  = zeros(T,R);

for r = 1:R
    for s = 1:NS
        samples(:,s,r) = normrnd( allpredsm(:,s,r), sqrt(allpredsv(:,s,r)));
    end
end
for r = 1:R
    wh_subjs = ~ismember(bestmodel(:,r), [1, 3]);
    meany(:,r)    = nanmean( tcn(:,wh_subjs,r), 2);
    vary(:,r)     = var(  tcn(:,wh_subjs,r), [], 2);
    samp_ypreds_mean(:,r) = nanmean( samples(:, wh_subjs, r), 2);
    samp_ypreds_var(:,r)  = var( samples(:, wh_subjs, r), [], 2);
    
    ypreds_mean(:,r) = nanmean(  models{ 1 }.allpredsm(:, wh_subjs, r), 2);
    ypreds_var(:,r)  = nanmean(  models{ 1 }.allpredsv(:, wh_subjs, r), 2);
end

meany = nanmean( meany, 2 );
vary  = nanmean( vary, 2 );
samp_ypreds_mean = nanmean( samp_ypreds_mean, 2 );
samp_ypreds_var = nanmean( samp_ypreds_var, 2 );
ypreds_mean = nanmean( ypreds_mean, 2 );
ypreds_var = nanmean( ypreds_var, 2 );

labels = [labels 'Other'];

% Plot the means
axes(hs(1));
plot( 1:T , ypreds_mean - 2, 'b-' , 'LineWidth' , lw ); hold on;
plot( 1:T , samp_ypreds_mean - 1, 'b--' , 'LineWidth' , lw ); hold on;
plot( 1:T , meany - 2, 'r-' , 'LineWidth' , lw );
title( titlestr );
xlim( [ 0 T ] );
ylim( [-M+1 NC]);

% Plot the var
axes(hs(2));
plot( 1:T , ypreds_var - 2 - 1, 'b-' , 'LineWidth' , lw ); hold on;
plot( 1:T , samp_ypreds_var - 2 - 1, 'b--' , 'LineWidth' , lw ); hold on;
plot( 1:T , vary - 2 - 1, 'r-' , 'LineWidth' , lw );
xlim( [ 0 T ] );
ylim( [-M+1 NC]);

labels = flip(labels);

axes(hs(1));
set(gca, 'YTick', -2:(NC-1));
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