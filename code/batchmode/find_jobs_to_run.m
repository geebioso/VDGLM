function [subjs_run, subjs_to_run] = find_jobs_to_run(whsim, isHPC, dotest)

if isHPC
    results_directory = '/pub/ggaut/VDGLM/Results';
else
    results_directory = fullfile( getenv('HOME'), 'Dropbox', 'FMRI', 'Projects',...
        'varianceGLM', 'Results');
end

input_directory = fullfile( results_directory, 'batch_analyses', 'single_jobs');
if dotest
    input_directory = [input_directory '_test'];
end

fprintf('Input directory: %s\n', input_directory);

%% Get list of files

files = dir(input_directory);
NF = length(files);
filenames = cell(NF,1);

for f = 1:NF
    filenames{f} = files(f).name;
end

% get rid of hidden files 
for f = 1:NF
    filenames{f} = files(f).name; 
    if filenames{f}(1) == '.'
      filenames{f} = [];  
   end
end

ii = cellfun( @(x) ~isempty(x), filenames, 'UniformOutput', true);
filenames = filenames(ii); 

% get only files that are from the appropriate simulation 
ii = cellfun( @(x)  contains( x, sprintf('whs%d', whsim)) , filenames, 'UniformOutput', true); 
filenames = filenames( ii ); 

NF = length(filenames); 

%% Get lists of subjects 
subjs_run = [];

for f = 1:NF
    filename = filenames{f};
    parts = split(filename, '_');
    subj_num =  str2num(parts{4});
    subjs_run = [subjs_run; subj_num];
end

% subjects that haven't been run yet
subjs_to_run = find( ~ismember( 1:875, subjs_run));


for s = 1:length(subjs_to_run)
   fprintf('%d\n', subjs_to_run(s));  
end
