

s=[];
for i = -1.8:0.2:0.2
   for j = 0.2:0.2:1.8
   
       [peak1,peak2,utotal,sink1,sink2] = moving_sink_inverse_model2_new(i,j);
       s1.mat = [peak1,peak2,utotal,sink1,sink2];
       s = [s,s1];
       save('two_inverse_sinks_mass.mat','s');
   end
   
end