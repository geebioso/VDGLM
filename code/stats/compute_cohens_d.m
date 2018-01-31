
whsims = [26];
input_directory = 'batchmode/Results_all';
whmodel = 1;
var_correction = 0;
dobraincontingency = 1;
doscatter = 1;
docompute = 1;
logscale = 0;
doabs = 0;

isHPC = 0;
dotest = 0;

addpath('ROI2NIfTI');
addpath('/Users/Garren/Dropbox/FMRI/BrainVisualization/NIfTI/NIfTI');
addpath('/Users/Garren/Dropbox/FMRI/Projects/varianceGLM/ROI2NIfTI/dicm2nii');

LOG = log4m.getLogger('');
LOG.setCommandWindowLevel(LOG.INFO);
LOG.setLogLevel(LOG.OFF);

[results_directory] = set_results_directory( isHPC );

%% Define Contrasts
contrasts = struct();
% contrasts.sim7.model1 = {[ 0 1 0 -1 0 0 0], [ 0 0 0 0 0 0 1], [ 0 0 0 0 0 1 0]};
% contrasts.sim7.model2 = {[ 0 1 0 -1 0 ]};
%
% contrasts.sim1.model1 = {[ 0 1 -1 0 0 0], [ 0 0 0 0 0 1]};
% contrasts.sim1.model2 = {[ 0 1 -1 0 0]};

contrasts.sim7.model1 = {
    [ -1 1 0 0 0 0 0], ...
    [ 0 0 0 0 -1 1 0], ...
    [ 0 1 0 0 0 0 0], ...
    [ 0 0 0 0 0 1 0], ...
    [ 1 0 0 0 0 0 0], ...
    [ 0 0 0 0 1 0 0], ...
    [ 0 0 1 0 0 0 0], ...
    [ 0 0 0 0 0 0 1]
    };

contrasts.sim7.names = {
    '2bck_minus_0bck_mean', ...
    '2bck_minus_0bck_var',...
    '2bck_minus_baseline_mean',...
    '2bck_minus_baseline_var', ...
    '0bck_minus_baseline_mean', ...
    '0bck_minus_baseline_var', ...
    'Instruction_minus_baseline_mean',...
    'Instruction_minus_baseline_var'
    };

contrasts.sim1.model1 = {[ 0 1 -1 0 0 0], [ 0 0 0 0 0 1]};
contrasts.sim1.names= {'2bck_minus_0bck_mean', '2bck_minus_baseline_variance'};

% same contrasts for pre-whitened and CIfTI analyses of the HCP data
contrasts.sim26 = contrasts.sim7;
contrasts.sim27 = contrasts.sim7;
contrasts.sim6 = contrasts.sim7;

D = struct();

%% Compute Cohen's D
if docompute
    for i  = 1:length(whsims)
        
        whsim = whsims(i);
        simfield = sprintf('sim%d', whsim);
        NC = length(fieldnames(contrasts.(simfield)));
        D.(simfield) = struct();
        D.(simfield).cohensd = cell(NC,1);
        D.(simfield).contrast_names = cell(NC,1);
        fprintf('Computing tests: sim %d, model %d\n', whsim, whmodel);
        
        %% Load Data
        
%         [models, allbicm, allllsm, bestmodelBIC, bestmodelCV, all_subjs, bad_subjs] = ...
%             load_results(results_directory, whsim, dotest, LOG);
        allparams = models{whmodel}.allparams;
        
        [NS, P, R] = size(allparams);
        
        nanp_subjs = any(any(isnan( allparams(:,:,:)), 3),2);
        if any( nanp_subjs)
            allparams(nanp_subjs,:,:) = [];
            %XS(nanp_subjs) = [];
            %FI(nanp_subjs,:) = [];
            I%FI(nanp_subjs,:) = [];
            %RCOND(nanp_subjs,:) = [];
            %tcn(:,nanp_subjs,:) = [];
            NS = NS - 1;
        end
        
        modelfield = sprintf('model%d', whmodel);
        NC = length(contrasts.(simfield).(modelfield)); % number of contrasts
        
        %% Cohen's D
        
        for c = 1:NC
            
            % specify contrast
            C = contrasts.(simfield).(modelfield){c};
            contrast_str = contrasts.(simfield).names{c};
            contrast_str = strrep(contrast_str,' ', '');
            g=sprintf('%d ', C);
            fprintf('\tContrast: %s, C = %s\n', contrast_str, g);
            
            % compute individual contrasts
            C2 = zeros(1,1,P);
            C2(:) = C;
            copes = sum( bsxfun( @times, permute(allparams, [1 3 2]), C2), 3); % checked that this is the same result as using a for loop
            
            % group level contrasts
            group_copes = mean( copes, 1);
            
            % Cohen's d
            cohensd = group_copes./std(copes,[],1);
            
            % correct variance for 1200 subjects
            if var_correction
                cohensd = sqrt( NS/1200 )*cohensd;
            end
            
            % Save to Structure
            D.(simfield).cohensd{c} = cohensd;
            D.(simfield).contrast_names{c} = contrast_str;
            % create nii file
            
        end
    end
    filename = fullfile(results_directory, '../ROI2NIfTI', 'files', ...
        sprintf('cohensd_whs%d_whmodel%d', whsim, whmodel));
    
    % output contrast names 
    f = fopen( sprintf('contrast_names_whs%d.txt', whsim), 'w'); 
    for c = 1:NC
        fprintf(f, '%s\n', D.(simfield).contrast_names{c});  
    end

    if whsim < 26
        error('YOU NEED TO WRITE CODE FOR COHENS D FOR THIS SIMULATION');
    else
        ROI2dscalar_nii_multi(whsim, D.(simfield).cohensd, ...
            D.(simfield).contrast_names, filename, 'pct');
    end
end

% options = struct();
% filename = fullfile(results_directory, '../ROI2NIfTI', 'files', sprintf('%s_contrast_cohensd_whs%d_whmodel%d.nii.gz', ...
%     contrast_str , whsim, whmodel));
% if whsim < 26
%     viewer = 'nii_viewer';
%     h = ROI2NIfTI(cohensd, filename, viewer, options);
% else
%     ROI2dscalar_nii(cohensd, filename, 'pct');
% end

%% Compare Mean and Variance Affects
if dobraincontingency
    
    for j = 1:length(whsims)
        
        whsim = whsims(j);
        BC = struct();
        cont_strs = {'2bck_minus_0bck', '2bck_minus_baseline', '0bck_minus_baseline'}; % contrast strings
        
        NC = length(cont_strs);
        f = fopen( sprintf('effect_names_whs%d.txt', whsim), 'w'); 
        for c = 1:NC
            fprintf(f, '%s\n', cont_strs{c}); 
        end
        effect_thresholds = [0.2, 0.5, 0.8];
        simfield = sprintf('sim%d', whsim);
        BC.(simfield) = struct(); 
        
        for i = 1:length(effect_thresholds)
            
            effect_threshold = effect_thresholds(i);
            
            switch effect_threshold
                case 0.2
                    effect_str = 'small';
                case 0.5
                    effect_str = 'medium';
                case 0.8
                    effect_str = 'large';
            end
            
            BC.(simfield).(effect_str) = cell(NC,1); 
            BC.(simfield).contrast_names = cont_strs; 
            
            for c = 1:NC
                
                fprintf('Processing %s at effect size %2.2f\n', cont_str, effect_threshold);
                
                % get contrast names
                cont_str = cont_strs{c};
                mean_cont_str = [cont_str '_mean'];
                var_cont_str = [cont_str '_var'];
                
                % get cohens d for the mean and variance
                mean_idx = strcmp( D.(simfield).contrast_names, mean_cont_str);
                var_idx = strcmp( D.(simfield).contrast_names, var_cont_str);
                mean_cohens = D.(simfield).cohensd{mean_idx};
                var_cohens = D.(simfield).cohensd{var_idx};
                
                % compute which effects occur
                mean_effect = and( abs(mean_cohens) > effect_threshold, abs(var_cohens)  < effect_threshold);
                var_effect= and( abs( var_cohens) > effect_threshold, abs(mean_cohens) < effect_threshold);
                both_effect= and( abs(mean_cohens) > effect_threshold, abs(var_cohens)  > effect_threshold);
                neither_effect = and( abs(mean_cohens) < effect_threshold, abs(var_cohens)  < effect_threshold );
                
                % store effects as categorical
                all_effects = zeros(R,1);
                all_effects( mean_effect )    = 1;
                all_effects( var_effect )     = 2;
                all_effects( both_effect )    = 3;
                all_effects( neither_effect ) = 0;
                
                BC.(simfield).(effect_str){c} = all_effects; 
            end
            
            filename = fullfile(results_directory, '../ROI2NIfTI', 'files', ...
                sprintf('effect_cohensd_%s_whs%d_whmodel%d', effect_str, whsim, whmodel));
            
            if whsim < 26
                error('YOU NEED TO WRITE CODE FOR COHENS D FOR THIS SIMULATION');
            else
                ROI2dscalar_nii_multi(whsim, BC.(simfield).(effect_str), ...
                    BC.(simfield).contrast_names, filename, 'pct');
            end
            
        end
    end
end

for i = {'small', 'medium', 'large'}
    fprintf('%s\n', i{1}); 
    for c = 1:NC
       fprintf('\tlength unique = %d\n', length(unique(BC.sim26.(i{1}){c})));  
    end
end
%% Mean Variance Scatter Plot

if doscatter
    cont_strs = { '2bck_minus_baseline', '0bck_minus_baseline', '2bck_minus_0bck' }; % contrast strings
    NC = length(cont_strs);
    effect_thresholds = [0.2, 0.5, 0.8];
    
    %% Plot
    
    f = figure(1); clf
    f.Position = [-7 387 1261 318];
    
    for c = 1:NC
        cont_str = cont_strs{c};
        mean_cont_str = [cont_str '_mean'];
        var_cont_str = [cont_str '_var'];
        mean_idx = strcmp( D.(simfield).contrast_names, mean_cont_str);
        var_idx = strcmp( D.(simfield).contrast_names, var_cont_str);
        mean_rois = D.(simfield).cohensd{mean_idx};
        var_rois = D.(simfield).cohensd{var_idx};
        
        cont_str = strrep(cont_str, '2', 'two');
        cont_str = strrep(cont_str, '0', 'zero');
        
        % effect sizes
        thresholds = [0.2, 0.5, 0.8];
        linetypes = {'-k', '--k', ':k' };
        
        groups = zeros(size(mean_rois));
        map = brewermap(4, 'Set1');
        for i = 1:length(thresholds)
            threshold = thresholds(i);
            ii = or( abs(mean_rois) > threshold, abs(var_rois) > threshold) ;
            groups(ii) = i;
        end
        
        % add negative effecdt sizes
        thresholds = [-thresholds, thresholds];
        linetypes = [linetypes, linetypes];
        
        h(c) = subplot(1,NC,c);
        
        if doabs
            x = abs(mean_rois);
            y = abs(var_rois);
        else
            x = mean_rois;
            y = var_rois;
        end
        
        % sc = scatter( x, y ); hold on;
        sc = gscatter( x, y, groups, map); hold on;
        
        if logscale
            h(c).XScale = 'log';
            h(c).YScale = 'log';
        else
            h(c).XScale = 'linear';
            h(c).YScale = 'linear';
        end
        
        h(c).FontSize = 13;
        
        xlabel('Mean d');
        ylabel('Var d');
        tit = strrep(cont_str, 'two', '2-');
        tit = strrep(tit, 'zero', '0-') ;
        tit = strrep(tit, 'bck', 'back');
        tit = strrep(tit, 'minus', '-');
        tit = strrep(tit, 'baseline', 'Fixation');
        title(strrep(tit, '_', ' '));
        
        
        % sc(1).Parent.Legend.Visible = 'on';
        
        for i = 1:length(thresholds)
            linetypenow = linetypes{i};
            % plot( [0.001 10], [thresholds(i), thresholds(i)], linetypenow); hold on;
            li = plot( h(c).XLim, [thresholds(i), thresholds(i)], linetypenow, 'LineWidth', 0.1 ); hold on;
        end
        
        for i = 1:length(thresholds)
            linetypenow = linetypes{i};
            % plot( [thresholds(i), thresholds(i)], [0.001 1], linetypenow); hold on;
            li = plot( [thresholds(i), thresholds(i)], h(c).YLim, linetypenow, 'LineWidth', 0.1 ); hold on;
        end
        
        legend(sc, {'No', 'Small', 'Medium', 'Large'});
        uistack(sc, 'top');
        
    end
    
    filename = fullfile( results_directory, '../images', sprintf('cohensd_mean_var_scatter_whs%d', whsim));
    print(filename, '-dpng');
    
    %% Plot Num regions by threshold
    h = figure(2); clf;
    h.Position = [-7 387 1261 318];
    
    for c = 1:NC
        cont_str = cont_strs{c};
        cont_str = strrep(cont_str, '2', 'two');
        cont_str = strrep(cont_str, '0', 'zero');
        
        mean_rois = abs(out.(cont_str).mean_rois);
        var_rois = abs(out.(cont_str).var_rois);
        
        mean_rois = sort(mean_rois);
        var_rois = sort(var_rois);
        
        thresholds = sort( [mean_rois; var_rois]);
        NT = length(thresholds);
        
        pct_mean = zeros(NT,1);
        pct_var = zeros(NT,1);
        pct_both = zeros(NT,1);
        for t = 1:NT
            
            mean_over = mean_rois > thresholds(t);
            var_over  = var_rois > thresholds(t);
            both_over = and( mean_rois > thresholds(t) , var_rois > thresholds(t) );
            mean_over(both_over) = false;
            var_over(both_over) = false;
            
            pct_mean(t) = sum( mean_over )/NT;
            pct_var(t)  = sum( var_over  )/NT;
            pct_both(t) = sum( both_over )/NT;
        end
        
        subplot(1,NC, c);
        plot(thresholds, pct_mean*100); hold on;
        plot(thresholds, pct_var*100 ); hold on;
        plot(thresholds, pct_both*100); hold on;
        xlabel('d');
        ylabel('pct rois');
        title(strrep(cont_str, '_', '\_'));
        
        legend('Location', 'best', {'Mean', 'Var', 'Both'});
    end
    
    filename = fullfile( results_directory, '../images', sprintf('cohend_num_above_threshold_whs%d', whsim));
    print(filename, '-dpng');
end

