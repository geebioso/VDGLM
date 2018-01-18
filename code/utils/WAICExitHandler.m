function generalExitHandler(src,event,filename)
   fprintf('\n');
   fprintf('Listener notified!\n');
   fprintf('Stan finished. Chains exited with exitValue = \n');
   disp(src.exit_value)
   fprintf('\n');
   
   samples = src.sim.samples; 
   
   
end