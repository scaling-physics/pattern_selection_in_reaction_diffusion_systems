filename = '../data/raw_data/SMCmodel3_sweep_c1.mat'; %datafile 
% System Parameters
L=2;
m = 0;
x = 0:0.01:L;

k1=1.5e-4;
k_minus1=3.6;
C=300; %constant total concentration
koff_arr = 0.001:0.005:0.15;
l = length(koff_arr);
%% 
nmax = 202; %loop over

j=1;
mode = zeros(l,nmax);
amp = zeros(l,nmax);
centroid_rate = zeros(l,nmax);
for k = 1:1:l  %flow control!
x = 1/koff_arr(k);
t = 0:1:60*60*x;
koff = koff_arr(k);
a=k1/k_minus1*C^2;
b=(koff+k_minus1)/k_minus1;
mode(j,1)=b;
amp(j,1) = b;
centroid_rate(j,1) = b;
mode(j,2) = a;
amp(j,2) = a;
centroid_rate(j,2) = a;
parfor i=3:nmax
[u1,u2] = SMCmodel3_out(C,koff,t,x);
%Finding fastest growing mode
for n=0:40
an(n+1) = 2*trapz(x,u2(end,:).*cos(pi*n*x/L)/1000)/L;
end
an = an(2:16);
[val,idx] = max(abs(an));
mode(j,i) = idx;
amp(j,i) = an(idx);

%Centroid
if idx == 2
centroid = diag(1./sum(u2,2))*u2*x';
I = get_last_monotonic(gradient(centroid),7);
ft = fittype('a+(c-a)*exp(-b*t)', 'independent', 't', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [L/2-0.5 0 0];
opts.StartPoint = [1 koff centroid(I)];
opts.Upper = [L/2+0.5 Inf L];

% Fit model to data.
[fitresult, gof] = fit( t(I:end)'-t(I), centroid(I:end), ft, opts );
rate = coeffvalues(fitresult);
rate = rate(2);
centroid_rate(j,i) = rate;
end
end

j=j+1;
end
save(filename,'nmax','L','C','x','t','mode', 'amp', 'centroid_rate') %saving data


