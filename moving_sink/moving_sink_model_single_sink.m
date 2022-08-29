function moving_sink_model_single_sink()

n_sinks=1;
D=0.3;
c=300;
mu=1;
delta=log(2)/50;
lambda = mu/delta;
L=1;
kappa=L*sqrt(delta/D)
v1=0.01;
v2=0.01;

Bp=[];
dG=[];
Gp=[];


G1=@(x,xi) kappa/2*(cosh(kappa*(x+xi)/L) + cosh(kappa*(abs(x-xi)-L)/L))/sinh(kappa);
Gx=@(x,xi) kappa^2/(2*L)*(sinh(kappa*(x+xi)/L) + sinh(kappa*(abs(x-xi)-L)/L).*sign(x-xi))/sinh(kappa) ;
Gx1 =@(x1) 2*kappa^2/(2*L)*(sinh(kappa*(2*x1)/L))/sinh(kappa);

t=0:0.1:60*10;

yinit=[0.3];


options=odeset('MaxStep',0.1);

[t,y1]=ode15s(@odes,t,yinit,options);
[t,y2]=ode15s(@odes1,t,yinit,options);

figure(1)
plot(t,y1(:,1),t,y2(:,1))
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
     

     G(1,1) = G1(y(1),y(1));   
     M=lambda*G+eye(n_sinks);
     mup=linsolve(M,lambda*c*ones(n_sinks,1));
     
     N = lambda*Gx1(y(1))*mup;
     mupd = linsolve(M,N);
     dydt = -v1*0.5*D*mupd;     
     
    end


 function dydt1=odes1(~,y)
     

     G(1,1) = G1(y(1),y(1));   
     M=lambda*G+eye(n_sinks);
     mup1=linsolve(M,lambda*c*ones(n_sinks,1));
     dG(1,1)=Gx(y(1)+0.00001,y(1));
     
     dydt1 = -v2*D*(dG+kappa^2*eye(n_sinks)/(2*L))*mup1; 

  
     
    end


end