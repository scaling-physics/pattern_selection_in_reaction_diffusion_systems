function moving_sink_model_two_sinks()

n_sinks=2;
D=0.3;
c=1;
mu=100;
delta=log(2)/50;
lambda = mu/delta;
L=1;
kappa=L*sqrt(delta/D)
v1=2;
v2=1;




G1 = @(x,xi) kappa/2*(cosh(kappa*(x+xi)/L) + cosh(kappa*(abs(x-xi)-L)/L))/sinh(kappa);
Gx = @(x,xi) kappa^2/(2*L)*(sinh(kappa*(x+xi)/L) + sinh(kappa*(abs(x-xi)-L)/L).*sign(x-xi))/sinh(kappa) ;

%derivative of G(x1,x1) w.r.t. x1
Gx1  =@(x1)     2*kappa^2/(2*L)*(sinh(kappa*(2*x1)/L))/sinh(kappa);

%derivative of G(x1,x2) w.r.t. x1
Gx12 =@(x1,x2)  kappa^2/(2*L)*(sinh(kappa*(x1+x2)/L)+sinh(kappa*(abs(x1-x2)-L)/L).*sign(x1-x2))/sinh(kappa);

%derivative of G(x2,x1) w.r.t. x1
Gx21 = @(x2,x1) kappa^2/(2*L)*(sinh(kappa*(x2+x1)/L)-sinh(kappa*(abs(x2-x1)-L)/L).*sign(x2-x1))/sinh(kappa);

 
t=0:0.1:1000;

yinit=[-0.3,0.4];

options=odeset('MaxStep',0.1);

[t,y1]=ode15s(@odes,t,yinit,options);
[t,y2]=ode15s(@odes1,t,yinit,options);
% 
figure(1)
plot(t,y1(:,1),t,y1(:,2),t,y2(:,1),t,y2(:,2))
ylim([-L/2,L/2])

% x=-L/2:0.01:L/2;
% for i=1:length(t)
%     for j=1:n_sinks
%         Gp(j,:)=1-G(y(i,j),y(i,:));
%     end
%     
%     M=b/delta*Gp+eye(n_sinks);
%     bp=linsolve(M,b*c*ones(n_sinks,1));
%     Bp=[Bp bp ];
%     
%     A(i,:)=c-(bp(1)+bp(2))/delta+G(x',y(i,1))*bp(1)/delta+G(x',y(i,2))*bp(2)/delta;
% end
% 
% figure(2)
% imagesc(t,x,A')
% colorbar;
% set(gca,'YDir','normal');
% 
% figure(3)
% plot(t,Bp(1,:),t,Bp(2,:))


    function dydt=odes(~,y)
     
     dydt = zeros(2,1);
             
     for j=1:n_sinks 
     G(j,:)=G1(y(j),y);
     end
     
     M=lambda*G+eye(n_sinks);
     mup=linsolve(M,lambda*c*ones(n_sinks,1));
     
     Gxm1(1,1) = Gx1(y(1));
     Gxm1(1,2) = Gx12(y(1),y(2));
     Gxm1(2,1) = Gx21(y(2),y(1));
     Gxm1(2,2)=  0;
     
     Gxm2(1,1) = 0;
     Gxm2(1,2) = Gx21(y(1),y(2));
     Gxm2(2,1) = Gx12(y(2),y(1));
     Gxm2(2,2) = Gx1(y(2));
     
     N1 = -lambda*Gxm1*mup;
     N2 = -lambda*Gxm2*mup;
     mupd1 = sum(linsolve(M,N1));
     mupd2 = sum(linsolve(M,N2));
     dydt(1) = 0.5*v1*D*mupd1;%-d/dx1 M=d/x1 sum mu_1
     dydt(2) = 0.5*v1*D*mupd2;
     
    end


    
     function dydt1=odes1(~,y)
     

         for j=1:n_sinks
            G(j,:)=G1(y(j),y);

            dG(j,:)=Gx(y(j)+0.00001,y);
         end
         M=lambda*G+eye(n_sinks);
         mup=linsolve(M,lambda*c*ones(n_sinks,1));
         dydt1 = -v2*D*(dG+kappa^2*eye(n_sinks)/(2*L))*mup; 
     end







end