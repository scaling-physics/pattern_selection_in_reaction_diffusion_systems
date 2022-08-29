%% Parameters
m = 0;
L= 1;
x = 0:0.01:L;
t = 0:1:60*60;
Du = 1;
Dv = 0.01;
% delta = 0.5;
% k_minus1=1;


% a = 1;
Gamma = 900;
d = Du/Dv;
% b = delta/k_minus1;


%%
mat = zeros(1,5);

figure
v = 0:1:100;  % plotting range from -5 to 5
[b a] = meshgrid(v); 
cond1=a+b>1;
d = 100;
cond2 = d.*(1-b) -a-2*sqrt(d.*a.*b) >0 ;

cond1 = double(cond1);  % convert to double for plotting
cond2 = double(cond2);
cond1(cond1 == 0) = NaN;  % set the 0s to NaN so they are not plotted
cond2(cond2 == 0) = NaN;
cond = cond1.*cond2;  % multiply the two condaces to keep only the common points
surf(b,a,cond)
view(0,90)


v = 0:0.1:5;
[x,y] = meshgrid(v);  % create a grid
ineq = x + y >= 1;    % some inequality
f = double(ineq);
figure
surf(x,y,f);
view(0,90) 

ezplot(@(b,a) cond1(b,a));
hold on
h = ezplot(@(b,a) cond2(a,b,d));


for m = 1:1:200
    for j = 1:1:50
        a=j*1;
        kappa= 0.1*m;
        b = kappa^2*2/Gamma;
%         kappa = sqrt(Gamma/(2*d));
        u0 = 1/a;
        v0 = 1;

        if cond(a,b) > 0 && cond1(a,b,d) > 0
             
      
            rr1=randn(size(x));
            rr2=randn(size(x));
            rr2=sum(rr1)*rr1/sum(rr2);
            IC = [u0*(1+0.1*rr1);v0*(1-0.1*rr2)];


            % Jacobian
            fu = -a;
            fv = -1;
            gu = a;
            gv = 1-b;



            nrange = [0:100*L];
            for i=1:length(nrange)
                %k=pi*n/L;
                k=pi*nrange(i);
                k_2 = k.^2;
                AA=Gamma.*[fu, fv;
                gu, gv] - [d*k_2 0; 0 k_2];
                lambda=eig(AA);
                [RE,I]=max(real(lambda));
                IM=imag(eig(I));
                A(i) = RE;
                A_im(i) = IM; 
            end


            
            nmodes = numel(A(A>1e-138));
            A(A<0) = 0;
            nmax = 0;
            
            if max(A) > 1e-138
            nmax = find(A==max(A))-1;
            end

            mat = [mat; kappa a nmodes nmax b];  

%             mat = [mat; Gamma a nmodes nmax];   
        
        end
       
    end
end

filename = "../data/raw_data/brusselator_turing_space_overlay_ldr_1.mat";
save(filename);