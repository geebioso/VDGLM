
whsims = [26];
input_directory = 'batchmode/Results_all';
whmodel = 1;
var_correction = 0;
dobraincontingency = 0; 
doscatter = 1;
docompute = 1;
logscale = 0;
doabs = 0;

isHPC = 0;
dotest = 0;

if whmodel ==2
   dobraincontingency = 0;  
   doscatter = 0; 
end

LOG = log4m.getLogger('test_log.txt');
LOG.setCommandWindowLevel(LOG.INFO);
LOG.setLogLevel(LOG.OFF);

[results_directory] = set_results_directory( isHPC, set_up_directory_structure );

%% Define Contrasts
contrasts = struct();
contrasts.sim7.model2.vec = {
    [ -1 1 0 0 ], ...    
    [  0 1 0 0 ],...
    [  1 0 0 0 ]
    };
contrasts.sim7.model2.names = {
    '2bck_minus_0bck_mean', ... 
    '2bck_minus_baseline_mean',...
    '0bck_minus_baseline_mean',
    }; 

contrasts.sim7.model1.vec = {
    [ -1 1 0 0 0 0 0], ...
    [ 0 0 0 0 -1 1 0], ...
    [ 0 1 0 0 0 0 0], ...
    [ 0 0 0 0 0 1 0], ...
    [ 1 0 0 0 0 0 0], ...
    [ 0 0 0 0 1 0 0], ...
    [ 0 0 1 0 0 0 0], ...
    [ 0 0 0 0 0 0 1]
    };

contrasts.sim7.model1.names = {
    '2bck_minus_0bck_mean', ...
    '2bck_minus_0bck_var',...
    '2bck_minus_baseline_mean',...
    '2bck_minus_baseline_var', ...
    '0bck_minus_baseline_mean', ...
    '0bck_minus_baseline_var', ...
    'Instruction_minus_baseline_mean',...
    'Instruction_minus_baseline_var'
    };

contrasts.sim1.model1.vec   = {[ 0 1 -1 0 0 0], [ 0 0 0 0 0 1]};
contrasts.sim1.model1.names = {'2bck_minus_0bck_mean', '2bck_minus_baseline_variance'};

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
        
        [models, allbicm, allllsm, bestmodelBIC, bestmodelCV, all_subjs, bad_subjs] = ...
            load_results(results_directory, whsim, dotest, LOG, 'analyze');
         
        allparams = models{whmodel}.allparams;
        
        [NS, P, R] = size(allparams);
        
        nanp_subjs = any(any(isnan( allparams(:,:,:)), 3),2);
        if any( nanp_subjs)
            allparams(nanp_subjs,:,:) = [];
            %XS(nanp_subjs) = [];
            %FI(nanp_subjs,:) = [];
            %FI(nanp_subjs,:) = [];
            %RCOND(nanp_subjs,:) = [];
            %tcn(:,nanp_subjs,:) = [];
            NS = NS - 1;
        end
        
        modelfield = sprintf('model%d', whmodel);
        NC = length(contrasts.(simfield).(modelfield).vec); % number of contrasts
        
        %% Cohen's D
        A = permute(allparams, [1 3 2]);
        for c = 1:NC
            
            % specify contrast
            C = contrasts.(simfield).(modelfield).vec{c};
            contrast_str = contrasts.(simfield).(modelfield).names{c};
            contrast_str = strrep(contrast_str,' ', '');
            g=sprintf('%d ', C);
            fprintf('\tContrast: %s, C = %s\n', contrast_str, g);
            
            % compute individual contrasts
            C2 = zeros(1,1,P);
            C2(:) = C;
            copes = sum( bsxfun( @times, A, C2), 3); % checked that this is the same result as using a for loop
            
            % group level contrasts
            group_copes = mean( copes, 1);
            
            % Cohen's d
            cohensd = group_copes./std(copes,[],1);
            
            % Non-paired Cohen's D 
%             s1 = std(A(:,:,1)); 
%             s2 = std(A(:,:,2)); 
%             T = 395; 
%             cohensd = group_copes./sqrt( (s1.^2 + s2.^2)*( T/(T-9) )); 
%             
          
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
    
    % save D
    filename = fullfile( results_directory, 'single_analyses', ...
        sprintf('cohensd_whs%d_whmodel%d', whsim, whmodel));
    
    save( filename, 'D');
    
    filename = fullfile(results_directory, '../ROI2NIfTI', 'files', ...
        sprintf('cohensd_whs%d_whmodel%d', whsim, whmodel));
    
    % output contrast names
    f = fopen( fullfile('..', 'plotting', sprintf('contrast_names_whs%d_mdl%d.txt', whsim, whmodel)), 'w');
    for c = 1:NC
        fprintf(f, '%s\n', D.(simfield).contrast_names{c});
    end
    
    if whsim < 26
        error('YOU NEED TO WRITE CODE FOR COHENS D FOR THIS SIMULATION');
    else
        
        temp = zeros(R, NC); 
        for c = 1:NC
           temp(:,c) = D.(simfield).cohensd{c};  
        end
        
        ROI2dscalar_nii_v2(temp, filename, D.(simfield).contrast_names, 'pct', dotest);
        
        save( fullfile( results_directory, '..', 'code', 'stats',...
            sprintf('cohensd_whs%d_whmodel%d.mat', whsim, whmodel)), 'D'); 
        
        % save fixed min/max bounds for mean and variance contrats 
        ismean = cellfun( @(x) contains( x, 'mean'), D.sim26.contrast_names); 
        isvar = ~ismean; 
        is_inst = cellfun( @(x) contains( x, 'Instruction'), D.sim26.contrast_names); 
       
        p = 1;
        idx = and(ismean, ~is_inst);
        upper = quantile( [D.sim26.cohensd{idx}], 1 - ( 1 - p )/2); 
        lower = quantile( [D.sim26.cohensd{idx}], (1-p)/2 );
        filenm =  fullfile('..', 'plotting', sprintf('bounds_whs%d_mdl%d_mean.txt', whsim, whmodel)); 
        f = fopen(filenm, 'w'); 
        fprintf(f, '%2.2f\n%2.2f\n', lower, upper); 
        fclose(f); 
        
        if whmodel == 1
            idx = and(isvar, ~is_inst);
            upper = quantile( [D.sim26.cohensd{idx}], 1 - ( 1 - p )/2);
            lower = quantile( [D.sim26.cohensd{idx}], (1-p)/2 );
            filenm =  fullfile('..', 'plotting', sprintf('bounds_whs%d_mdl%d_var.txt', whsim, whmodel));
            f = fopen(filenm, 'w'); 
            fprintf(f, '%2.2f\n%2.2f\n', lower, upper); 
            fclose(f); 
        
        end
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
                % get contrast names
                cont_str = cont_strs{c};
                fprintf('Processing %s at effect size %2.2f\n', cont_str, effect_threshold);
                
                % get cohens d for the mean and variance
                mean_cont_str = [cont_str '_mean'];
                mean_idx = strcmp( D.(simfield).contrast_names, mean_cont_str);
                mean_cohens = D.(simfield).cohensd{mean_idx};
                
                var_cont_str = [cont_str '_var'];
                var_idx = strcmp( D.(simfield).contrast_names, var_cont_str);
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
                ROI2dscalar_nii_multi(BC.(simfield).(effect_str), ...
                    BC.(simfield).contrast_names, filename, 'pct');
            end
            
        end
    end
    
    
    for i = {'small', 'medium', 'large'}
        fprintf('%s\n', i{1});
        for c = 1:NC
            fprintf('\tlength unique = %d\n', length(unique(BC.sim26.(i{1}){c})));
        end
    end
    
end
%% Mean Variance Scatter Plot

if doscatter
    cont_strs = { '2bck_minus_baseline', '0bck_minus_baseline', '2bck_minus_0bck' }; % contrast strings
    NC = length(cont_strs);
    effect_thresholds = [0.2, 0.5, 0.8];
    
    %% Get min and max cohensd for mean and var 
    ii = cellfun( @(x) ~isempty(x), strfind( D.sim26.contrast_names, 'mean')); 
    MD = cell2mat(D.sim26.cohensd(ii)); 
    xlims = [min(MD(:)), max(MD(:))]; 
    
    VD = cell2mat(D.sim26.cohensd(~ii)); 
    ylims = [min(VD(:)), max(VD(:))]; 
    
    delta = 0.1; 
    
    %% Plot
    
    f = figure(1); clf
    f.Position = [-7 387 1261 318];
    
    for c = 1:NC
        cont_str = cont_strs{c};
        
        mean_cont_str = [cont_str '_mean'];
        mean_idx = strcmp( D.(simfield).contrast_names, mean_cont_str);
        mean_rois = D.(simfield).cohensd{mean_idx};
        mean_groups = zeros(size(mean_rois));
        
        var_cont_str = [cont_str '_var'];
        var_idx = strcmp( D.(simfield).contrast_names, var_cont_str);
        var_rois = D.(simfield).cohensd{var_idx};
        var_groups = zeros( size(var_rois)); 
        
        cont_str = strrep(cont_str, '2', 'two');
        cont_str = strrep(cont_str, '0', 'zero');
        
        % effect sizes
        thresholds = [0.2, 0.5, 0.8];
        linetypes = {'-k', '--k', ':k' };
        
        map = brewermap(9, 'Set1');
        temp = map(9,:); 
        map(2:9,:) = map(1:8,:);
        map(1,:) = temp;
        for i = 1:length(thresholds)
            if i == 1
                lower = 0;
            else
                lower = thresholds(i-1);  
            end
            upper = thresholds(i);
            
            if i == 1
               no_effect_idx = and( abs(mean_rois) < upper, abs(var_rois) < upper); 
            end
            ii = and( abs(mean_rois) > lower, abs(mean_rois) < upper) ;
            mean_groups(ii) = i;
            
            ii = and( abs(var_rois) > lower, abs(var_rois) < upper) ;
            var_groups(ii) = i;
        end
        
        % get which regions have mean and variance effect sizes in the same
        % threshold bin 
        groups = (var_groups == mean_groups) + 1; 
        groups(no_effect_idx) = 0; 
        
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
        
        h(c).XLim = [xlims(1) - delta, xlims(2) + delta];
        h(c).YLim = [ylims(1) - delta, ylims(2) + delta];
        
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
        
        legend(sc, {'no effect', 'diff size', 'same size'});
        uistack(sc, 'top');
        
    end
    
    filename = fullfile( results_directory, '../images', sprintf('cohensd_mean_var_scatter_whs%d', whsim));
    print(filename, '-dpng');
    
    %     %% Plot Num regions by threshold
    %     h = figure(2); clf;
    %     h.Position = [-7 387 1261 318];
    %
    %     for c = 1:NC
    %         cont_str = cont_strs{c};
    %         cont_str = strrep(cont_str, '2', 'two');
    %         cont_str = strrep(cont_str, '0', 'zero');
    %
    %         mean_rois = abs(out.(cont_str).mean_rois);
    %         var_rois = abs(out.(cont_str).var_rois);
    %
    %         mean_rois = sort(mean_rois);
    %         var_rois = sort(var_rois);
    %
    %         thresholds = sort( [mean_rois; var_rois]);
    %         NT = length(thresholds);
    %
    %         pct_mean = zeros(NT,1);
    %         pct_var = zeros(NT,1);
    %         pct_both = zeros(NT,1);
    %         for t = 1:NT
    %
    %             mean_over = mean_rois > thresholds(t);
    %             var_over  = var_rois > thresholds(t);
    %             both_over = and( mean_rois > thresholds(t) , var_rois > thresholds(t) );
    %             mean_over(both_over) = false;
    %             var_over(both_over) = false;
    %
    %             pct_mean(t) = sum( mean_over )/NT;
    %             pct_var(t)  = sum( var_over  )/NT;
    %             pct_both(t) = sum( both_over )/NT;
    %         end
    %
    %         subplot(1,NC, c);
    %         plot(thresholds, pct_mean*100); hold on;
    %         plot(thresholds, pct_var*100 ); hold on;
    %         plot(thresholds, pct_both*100); hold on;
    %         xlabel('d');
    %         ylabel('pct rois');
    %         title(strrep(cont_str, '_', '\_'));
    %
    %         legend('Location', 'best', {'Mean', 'Var', 'Both'});
    %     end
    %
    %     filename = fullfile( results_directory, '../images', sprintf('cohend_num_above_threshold_whs%d', whsim));
    %     print(filename, '-dpng');
end

