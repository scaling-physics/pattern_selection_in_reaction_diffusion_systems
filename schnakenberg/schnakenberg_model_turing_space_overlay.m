%% Parameters
m = 0;
L= 2;
x = 0:0.01:L;
t = 0:1:1000;


gamma = 150;
Du = 2;
Dv = 0.01;
d= Du/Dv;


%%
mat = zeros(1,4);

C = 1;
for j = 1:100
    for m = 1:100
        k2 = j;
        k3 = m;
        k4 = C*k2;
 
        Gamma = k2*L^2/Dv;
        a = k3/k2*C^2;
        u0 = C;
        v0 = k2/(k3*C);

        rr1 = randn(size(x));
        rr2 = randn(size(x));
        rr2 = sum(rr1)*rr1/sum(rr2);
        IC = [u0*(1+0.1*rr1);v0*(1-0.1*rr2)];
    

            % Jacobian
            fu=k2;
            fv=k3*C^2;
            gu=-2*k2;
            gv=-k3*C^2;



            nrange = [0:50*L];

            for i=1:length(nrange)
                %k=pi*n/L;
                k=pi*nrange(i)/L;
                k_2 = k.^2;
                AA=[fu-Du*k_2, fv;
                gu, gv-Dv*k_2];
                lambda=eig(AA);
                [RE,I]=max(real(lambda));
                IM=imag(eig(I));
                A(i) = RE;
                A_im(i) = IM; 
            end

    %         figure(8);
    %         plot(nrange,A,'-o',nrange,A_im)
    %         title('Dispersion Relation','FontSize',12);
    %         xlabel('n')
    %         ylabel('growth rate')
    %         set(gca,'FontName', 'Helvetica','FontSize', 10)



            nmodes = numel(A(A>0.0001));
            A(A<0) = 0;
            nmax = 0;

            if max(A) > 0.0001
            nmax = find(A==max(A))-1;
            end

%
%             mat = [mat; sqrt(Gamma/(2*d)) Gamma*a nmodes nmax];
            mat = [mat; sqrt(k2) k3 nmodes nmax];
%         end
    end
end

filename = "../data/raw_data/schnakenberg_turing_space_overlay_ldr.mat";
save(filename);