function [ Vpermute ] = permute_data( V, Nperm, S, K, R )
% Function to permute the data over time (keeping the same subject, region,
% and task)
% Nperm is number of permuted datasets to create 

T = size(V,4); 
Vpermute = zeros(Nperm, S, K, R, T); 


for n = 1:Nperm
    for s = 1:S
        for k = 1:K
            for r = 1:R
                Vpermute(n,s,k,r,:) = V(s,k,r, randperm(T)); 
            end
        end
    end
end

end

