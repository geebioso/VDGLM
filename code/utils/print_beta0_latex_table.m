function [ ] = print_beta0_latex_table( filename, samples, beta0)
% puts the confusion matrix into format for latex code

var_beta0 = squeeze(mean(mean(samples.var_beta0)));
std_beta0 = squeeze(mean(mean(samples.std_beta0)));

S = size(var_beta0, 1);
K = size(var_beta0, 2);


f = fopen(filename, 'w');

% headers 
fprintf(f, '\\begin{tabular} {c c r r r }\n'); 
fprintf(f, '{\\bf Subject } & {\\bf Task} & {\\bf True} & {\\bf VarGLM } & { \\bf StdGLM } \\\\ \\hline \\\\[-2ex]\n');  

for s = 1:S
    for k = 1:K
        fprintf(f, '%d & %d & %2.2f & %2.2f & %2.2f \\\\\n', s, k, beta0(s,k), var_beta0(s,k), std_beta0(s,k)); 
    
   
    end
end
fprintf(f, '\\\\ \\\\[-2ex] \\hline \\\\ \n\n');
fprintf(f, '\\end{tabular}'); 

fclose(f);

fprintf('\nVarGLM Beta0:\n'); 
disp(var_beta0)

fprintf('StdGLM Beta0:\n'); 
disp(std_beta0)

fprintf('True Beta0:\n'); 
disp(beta0); 




