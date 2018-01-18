function [subjs_run, subjs_to_run] = find_jobs_to_run(isHPC, dotest)

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

if NF > 2
    if strcmp(filenames{1}, '.')
        filenames = filenames(3:end);
        NF = NF - 2;
    end
end

if NF > 1
    if strcmp(filenames{1}, '.DS_Store')
        filenames(1) = [];
        NF = NF - 1;
    end
end

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
