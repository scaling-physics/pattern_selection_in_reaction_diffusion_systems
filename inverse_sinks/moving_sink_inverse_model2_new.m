function [peak1,peak2,sink1,sink2,utotal] = moving_sink_inverse_model2_new(x1,x2)
%solving as a DAE

n_sinks=2;

L=4;
beta=1.5e-4;
gamma=3.6;
delta=log(2)/50;
Du=0.3;
Dv=0.0012;
rho=6*delta*sqrt(Dv*(gamma+delta))/beta/L;
c=300;
v1=1;
kappa=L*sqrt(delta/Du)
sigma = rho/(delta*c^2);
dG=[];
Gp=[];


% G=@(x,xi) 1-kappa/2*(   cosh(kappa*(x+xi)/L)  +  cosh(kappa*(abs(x-xi)-L)/L)  )/sinh(kappa);
%Gx=@(x,xi) -kappa^2/(2*L)*(   sinh(kappa*(x+xi)/L)   +  sinh(kappa*(abs(x-xi)-L)/L).*sign(x-xi)   )/sinh(kappa) ;
G1 = @(x,xi) kappa/2*(cosh(kappa*(x+xi)/L) + cosh(kappa*(abs(x-xi)-L)/L))/sinh(kappa);
Gx = @(x,xi) kappa^2/(2*L)*(sinh(kappa*(x+xi)/L) + sinh(kappa*(abs(x-xi)-L)/L).*sign(x-xi))/sinh(kappa) ;

%derivative of G(x1,x1) w.r.t. x1
Gx1  =@(x1)     2*kappa^2/(2*L)*(sinh(kappa*(2*x1)/L))/sinh(kappa);

%derivative of G(x1,x2) w.r.t. x1
Gx12 =@(x1,x2)  kappa^2/(2*L)*(sinh(kappa*(x1+x2)/L)+sinh(kappa*(abs(x1-x2)-L)/L).*sign(x1-x2))/sinh(kappa);

%derivative of G(x2,x1) w.r.t. x1
Gx21 = @(x2,x1) kappa^2/(2*L)*(sinh(kappa*(x2+x1)/L)-sinh(kappa*(abs(x2-x1)-L)/L).*sign(x2-x1))/sinh(kappa);


t=0:0.1:60*100;
x1 = -1.8;
x2 = 1.4;

[rhoprimes,~,Atot ]=get_bprimes(rho,c,delta,Du,L,[x1,x2]);
rhoprimes.rhop1;
rhoprimes.rhop2;
rhoprimes.rhop1(imag(rhoprimes.rhop1)~=0)=0;
rhoprimes.rhop2(imag(rhoprimes.rhop1)~=0)=0;

[~,I2]=max(rhoprimes.rhop1+rhoprimes.rhop2);


%variables are x1, x2
yinit=[x1,x2,double(rhoprimes.rhop1(I2)),double(rhoprimes.rhop2(I2))];

xx=yinit(1:2);


M = [1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0];
options=odeset('Mass', M);

%%%get initial conditions for rho'

Atot;
%%%
% 
[t,y1]=ode15s(@odes,t,yinit,options);
[t,y]=ode15s(@odes1,t,yinit,options);




figure(1)
% clf;

% plot(t,y1(:,1:2))
plot(y1(:,1)/L,y1(:,2)/L)
hold on
plot(y1(:,1)/L,y1(:,2)/L)
% peak1 = y(:,1)/L;
% peak2 = y(:,2)/L;
% utotal = (c-(y1(:,3)+y1(:,4))/delta)/c;
% ylim([-1/2,1/2])
% xlim([-1/2,1/2])
% % 
% figure(2)
% clf;
% plot(t,utotal)
peak1=y(:,1)/L;
peak2= y(:,2)/L;
sink1 = y1(:,1)/L;
sink2 = y1(:,2)/L;
utotal = (c-(y1(:,3)+y1(:,4))/delta)/c;


    function dydt=odes(~,y)
        dydt=zeros(4,1);
        x=y(1:2);
        for j=1:n_sinks
            dG(j,:)=Gxp(x(j),x);
        end
        for j=1:length(x)
            Gp(j,:)=G1(x(j),x);
        end
%         [rhoprimes,~,Atot] = get_bprimes(rho,c,delta,Du,L,[y(1),y(2)]);
%         rhoprime = [max(real(rhoprimes.rhop1)); max(real(rhoprimes.rhop2))];
        dydt(1:2)= v1*Du/2*(2*dG*y(3:4)-y(3:4)*kappa^2/(L));
        dydt(3)=(c-Gp(1,1)*y(3)-Gp(1,2)*y(4))*y(3)-rho/delta;%eqn for rho'1
        dydt(4)=(c-Gp(2,1)*y(3)-Gp(2,2)*y(4))*y(4)-rho/delta;%eqn for rho'2

    end

    function out=Gxp(x,xi)%G(xj+;xi)
        out=size(xi);
        I=x<xi;
        out(I)=-kappa^2/(2*L)*(   sinh(kappa*(x+xi(I))/L)   -  sinh(kappa*(xi(I)-x-L)/L)   )/sinh(kappa) ;
        out(~I)=-kappa^2/(2*L)*(   sinh(kappa*(x+xi(~I))/L)   +  sinh(kappa*(x-xi(~I)-L)/L)   )/sinh(kappa) ;        
    end

     function dydt1=odes1(~,y)
     
     dydt1 = zeros(4,1);
     x=y(1:2);        
     for j=1:n_sinks 
     G(j,:)=G1(y(j),x);
     end
     
     
%      mup=linsolve(M,mu*c*ones(n_sinks,1));
%      [rhoprimes,~,Atot] = get_bprimes(rho,c,delta,Du,L,[y(1),y(2)]);
%      rhoprime = [max(real(rhoprimes.rhop1)); max(real(rhoprimes.rhop2))];
     
     rhoprime=y(3:4);
     M1 = -diag(rhoprime)*G; % row wise multiplication
     M2 = diag(c-G*rhoprime);
     M = M1+M2;
     Gxm1(1,1) = Gx1(y(1));
     Gxm1(1,2) = Gx12(y(1),y(2));
     Gxm1(2,1) = Gx21(y(2),y(1));
     Gxm1(2,2)=  0;
     
     Gxm2(1,1) = 0;
     Gxm2(1,2) = Gx21(y(1),y(2));
     Gxm2(2,1) = Gx12(y(2),y(1));
     Gxm2(2,2) = Gx1(y(2));
     
     
     N1 = -(Gxm1*rhoprime).*rhoprime;
     N2 = -(Gxm2*rhoprime).*rhoprime;
     rhopd1 = sum(linsolve(M,N1));
     rhopd2 = sum(linsolve(M,N2));

     dydt1(1) = n_sinks*0.5*v1*Du*rhopd1;%-d/dx1 M=d/x1 sum mu_1
     dydt1(2) = n_sinks*0.5*v1*Du*rhopd2;
     
     dydt1(3)=(c-G(1,1)*y(3)-G(1,2)*y(4))*y(3)-rho/delta;%eqn for rho'1
     dydt1(4)=(c-G(2,1)*y(3)-G(2,2)*y(4))*y(4)-rho/delta;%eqn for rho'2
     
    end



end