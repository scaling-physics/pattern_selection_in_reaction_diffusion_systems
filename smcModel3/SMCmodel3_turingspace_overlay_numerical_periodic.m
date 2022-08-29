filename = '../data/raw_data/SMCmodel3_sweep_turingspace_overlay_numerical_analysis_periodic_sweep.mat'; %datafile 
% Simulation settings
L=4;
m = 0;
x = 0:0.01:L;
% t= 1:1:60*60;

% Parameters
k_minus1=3.6;
C=300; %constant total concentration
d = 25;

cond=@(d,b,a) d*(b+1).*(a-b-1)./(a+b+1)+(-(a.^2+b.^2+4*a.*b+3.*a+b)./(a+b+1))-2*sqrt(d*(b).*(a+b+1));

%% 
nmax = 52; %loop over

mat = zeros(1,nmax);
mat_polarity = zeros(1,3);
j=1;
for k = 1:250  %flow control!
    for l =0:250

        b = 0.001*k;
        a = 0.1*l;              
       
        if cond(d,b,a) > 0
            koff = b*k_minus1;
            k1 = a*k_minus1/C^2;        
            mat(j,1)=b;
            mat(j,2)=a;
            

            parfor i=3:nmax
                
%             [u1,u2] = SMCmodel3_out(C,koff,k1,a,b,t,x);
            [u1,u2] = periodic_sweep(C,koff,k1,a,b,x);


            grad=diff(u2(end,:));
            grad = grad(abs(grad)>2);

            m = diff(sign(grad));
            if isempty(grad)
                npeaks = 0;
            else
                npeaks = length(find(m))/2+0.5;
            end
            mat(j,i) = npeaks;    
            
            end

            j=j+1;
        end
    end
end
% save(filename, 'L','C','x','t','mode', 'amp', 'centroid_rate') %saving data
save(filename);

