function moving_sink_model()

n_sinks=2;
D=0.01;
c=1;
b=1;
delta=0.002;
L=1;
kappa=L*sqrt(delta/D)
v=50;

Bp=[];
bp=[];
dG=[];
Gp=[];


G=@(x,xi) 1-kappa/2*(   cosh(kappa*(x+xi)/L)  +  cosh(kappa*(abs(x-xi)-L)/L)  )/sinh(kappa);
Gx=@(x,xi) -kappa^2/(2*L)*(sinh(kappa*(x+xi)/L) + sinh(kappa*(abs(x-xi)-L)/L).*sign(x-xi))/sinh(kappa) ;

t=0:0.1:1000;

yinit=[-0.3,0.4];


options=odeset('MaxStep',0.1);

[t,y]=ode15s(@odes,t,yinit,options);


figure(1)
plot(t,y)
ylim([-L/2,L/2])

x=-L/2:0.01:L/2;
for i=1:length(t)
    for j=1:n_sinks
        Gp(j,:)=1-G(y(i,j),y(i,:));
    end
    
    M=b/delta*Gp+eye(n_sinks);
    bp=linsolve(M,b*c*ones(n_sinks,1));
    Bp=[Bp bp ];
    
    A(i,:) = c-(bp(1)+bp(2))/delta+G(x',y(i,1))*bp(1)/delta+G(x',y(i,2))*bp(2)/delta;
end

figure(2)
imagesc(t,x,A')
colorbar;
set(gca,'YDir','normal');

figure(3)
plot(t,Bp(1,:),t,Bp(2,:))


    function dydt=odes(~,y)
        for j=1:n_sinks
            Gp(j,:)=1-G(y(j),y);
            dG(j,:)=Gx(y(j)+0.00001,y);
        end
        
     M=b/delta*Gp+eye(n_sinks);
     bp=linsolve(M,b*c*ones(n_sinks,1))
     
     dydt=v*(2*dG-kappa^2*eye(n_sinks)/(L))*bp; 
    end







end