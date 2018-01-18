function [ post_means ] = get_posterior_means( samples)

fieldnms = fieldnames(samples); 
NF = length(fieldnms); 
post_means = struct(); 

for i = 1:NF
   fieldnow = fieldnms{i};
   sampnow = samples.(fieldnow); 
   post_means.(fieldnow) = squeeze(mean(sampnow,1)); 
end


end

