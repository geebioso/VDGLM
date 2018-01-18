function [] = clear_samples(wd, filename, nchains)

split_fn =  split(filename, '.'); 
for i = 1:nchains
    
    fn = [wd split_fn{1} '-' num2str(i) split_fn{2}]; 
    delete(fn); 
   
end



end