filename = '/Users/srik/Documents/MATLAB/flow/data/raw_data/SMCmodel3_sweep_c300_6.mat'; %datafile 
load(filename);
data = table(mode,amp,centroid_rate);
range = length(koff_arr);
m=zeros(range,nmax);
mode_amp=zeros(range,nmax);

for i=1:range
mode_amp(i,1) = data.mode(i,1);
mode_amp(i,2) =  data.mode(i,2);
for j=3:nmax
if abs(data.amp(i,j)) > 1e-5
    if data.amp(i,j) < 0 
        mode_amp(i,j) = -data.mode(i,j);
    else
        mode_amp(i,j) =  data.mode(i,j);
    end
end
if data.amp(i,j)<0 & data.centroid_rate(i,j) <0.1
    m(i,j)=data.centroid_rate(i,j);
end
end
end

k_minus1=3.6;
edges = 0.5:1:8.5;
ma = zeros(range,6);
for i=1:range
N = histcounts(data.mode(i,3:end),edges);
[val,idx] = max(N(1,:));
ma(i,1) = data.mode(i,1);
ma(i,2) = data.mode(i,2);
ma(i,3) = k_minus1*(ma(i,1)-1);
ma(i,4) = 1/ma(i,3);
ma(i,5) = idx;
rate = sum(m(i,3:end)); 
len = length(find(m(i,3:end)));
rate = rate/len;
ma(i,6)=rate;
end
% colNames = {'b','a','koff','T','mode','rate'};
% mTable = array2table(ma,'VariableNames',colNames);



