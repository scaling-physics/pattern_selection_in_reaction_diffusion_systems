function [u1,u2]=SMCmodel3_out(C,koff,t,x)

k1 = 1.5e-4;
k_minus1 = 3.6;
kon = C*koff;

Du = 0.3;
Dv = 0.012;

a = k1/k_minus1*C^2;
b = koff/k_minus1;
v0 = C*a/(a+b+1);
u0 = C-v0;

rr1=randn(size(x));
rr2=randn(size(x));
rr2=sum(rr1)*rr1/sum(rr2);
IC = [u0*(1+0.1*rr1);v0-0.1*u0*rr2];



options=odeset('RelTol',1e-6,'AbsTol',1e-12);
sol = pdepe(0,@pdes,@ic,@bc,x,t,options);
u1 = sol(:,:,1);
u2 = sol(:,:,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------
function [c,f,s] = pdes(~,~,u,DuDx)

c = [1; 1];
f = [Du;Dv] .* DuDx; 
F = k1*u(1)*(u(1)+u(2))^2-k_minus1*u(2);
s = [-F+kon-koff*u(1);F-koff*u(2)]; 
end
% --------------------------------------------------------------
function out = ic(xx)

out = IC(:,x==xx);
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bc(~,~,~,~,~)
pl = [0;0]; %negative=flux out, positive= flux in
ql = [1; 1]; 
pr = [0;0]; % negative=flux in, positive=flux out
qr = [1; 1];
end


end