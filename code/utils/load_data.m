function [data, ND, S, K, R, V, T, modelfile] = load_data( data_filename, modelfile)

load(data_filename); % loads orignal data in variable samples 

data = expand_samples(samples); 

V = samples.V;
% V = permute(V,[2,3,4,1,5]); 


ND = size(V,1); 
S = size(V,2);                              % number of subject
K = size(V,3);                              % number of tasks
R = size(V,4);                              % number of regions
T = size(V,5); 

if R == 1
    modelfile = [modelfile '_one_region.stan']; 
    
else
    modelfile = [modelfile '.stan'];
end

V = squeeze(V);

end