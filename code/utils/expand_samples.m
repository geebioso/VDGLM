function [ new_samples ] = expand_samples( samples )
% Takes a structure array of samples and puts all samples ina structure

new_samples = struct();
fieldnms = fieldnames(samples); 
for i = 1:size(samples,2)
    for j = 1:length(fieldnms)
        fieldnow = fieldnms{j}; 
        temp = samples(i); 
        if i == 1 
            new_samples.(fieldnow) = temp.(fieldnow); 
        else
            new_samples.(fieldnow) = [new_samples.(fieldnow); temp.(fieldnow)]; 
        end
    end
end

end

