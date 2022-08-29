function [rhoprimes T Atot]=get_bprimes(rho,c,delta,Du,L,y)

kappa=L*sqrt(delta/Du);
G=@(x,xi) kappa/2*(   cosh(kappa*(x+xi)/L)  +  cosh(kappa*(abs(x-xi)-L)/L)  )/sinh(kappa);

for j=1:length(y)
     Gp(j,:)=G(y(j),y);
end


syms rhop1 rhop2

eqns=[rho==(c-Gp(1,1)*rhop1/delta-Gp(1,2)*rhop2/delta)*rhop1,rho==(c-Gp(2,1)*rhop1/delta-Gp(2,2)*rhop2/delta)*rhop2];
rhoprimes=vpasolve(eqns,[rhop1,rhop2]);

Atot=L*c-L*(rhoprimes.rhop1+rhoprimes.rhop2)/delta;
T=Atot(imag(Atot)==0);
if ~isempty(T)
    T=double(min(T));
else
    T=NaN;
end

end