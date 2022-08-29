%%%%%   in units of nM
%%%% u(u+v)^2
L= 2;
m = 0;
% x = 0:0.01:L;
% t = 0:1:60*60;

k1= 1.5e-4;
k_minus1=3.6;
Du=0.3;
Dv=0.012;
koff = log(2)/50;
C=300;
% kon=C*koff;


Gamma = k_minus1*L^2/Dv;
b=koff/k_minus1;
a=k1/k_minus1*C^2;
d = Du/Dv;

mat = zeros(1,4);

v0 = a/(a+b+1);
u0 = 1-v0;


fu = -a-2*a*u0-b;
fv = 1-2*a*u0;
gu = a+2*u0;
gv = 2*a*u0-1-b;



nrange = [0:4*L];
krange = 0:0.01:pi*nrange(end);

for i=1:length(nrange)
    %k=pi*n/L;
    k=pi*nrange(i);
    k2 = k.^2;
    AA=[Gamma*fu-d*k2, Gamma*fv;
        Gamma*gu, Gamma*gv-k2];
    lambda=eig(AA);
    [RE,I]=max(real(lambda));
    IM=imag(eig(I));
    A(i) = RE;
    A_im(i) = IM; 
end

%%alternative: explicit formula for the greater eigenvalue
%Eig=@(k2) ((fu+gv-k2*(Du+Dv))+sqrt((fu+gv-k2*(Du+Dv)).^2-4*(fu*gv-gu*fv-k2*(fu*Dv+gv*Du)+k2.^2*Du*Dv)   )  )/2;

figure(1);
plot(nrange,A,'-o',nrange,A_im)
title('Dispersion Relation','FontSize',12);
xlabel('n')
ylabel('growth rate')
set(gca,'FontName', 'Helvetica','FontSize', 10)




