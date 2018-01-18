
whsim = 7;
whmodel = 1;

%% Load AC and DW
load( fullfile( 'Results', sprintf('autocorrelation_whs%d_whm%d.mat', whsim, whmodel)));

%% Plot Mean AC with error bars

meanAC = squeeze(nanmean(nanmean( AC, 1), 2));

figure(1); clf;
plot(meanAC, 'bo', 'LineWidth', 1.5);
xlabel('lag');
ylabel('autocorrelation');

title('Mean Autocorrelation');

%% Plot AC at each Lag
figure(2); clf;

for i = 1:9
    subplot(3,3, i);
    imagesc(AC(:,:,i)); colorbar;
    title(sprintf('lag %d', i));
end

%% Plot DW at each lag

h = figure(3); clf;

imagesc(DW < 0.05); colorbar;
title('Durbin Watson Test For Autocorrelation');
ylabel('Subject'); 
xlabel('Region'); 

filename = fullfile('images', sprintf('dw_autocorr_test_whs%d_whm%d', whsim, whmodel)); 
print(filename, '-dpng'); 
