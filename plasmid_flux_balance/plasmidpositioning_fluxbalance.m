
L=1;


kp=1;%production rate
D=0.01;%diffusion constant
C=kp/(2*D);

figure(1)

xp=0.5; %plasmid position
xl=0:0.01:xp;
xr=xp:0.01:L;
subplot(3,1,1)
Al=C*(xp^2-xl.^2);
Ar=C*((L-xp)^2-(L-xr).^2);
plot(xl,Al,xr,Ar)
%ylim([0,1.1])


subplot(3,1,2)

dflux=[];
Atot=[];
for xp=0:0.01:1
xl=0:0.01:xp;
xr=xp:0.01:L;
Al=C*(xp^2-xl.^2);
Ar=C*((L-xp)^2-(L-xr).^2);
fluxL=D*abs(gradient(Al))/0.01;
fluxR=D*abs(gradient(Ar))/0.01;
dflux=[dflux fluxR(1)-fluxL(end)];
Atot=[Atot sum([Al,Ar])];
end

plot(0:0.01:1,dflux)
subplot(3,1,3)
plot(0:0.01:1,Atot)


figure(2)

%%with 2 sinks
xp1=0;
xp2=1;
x1=0:0.01:xp1;
x2=xp1:0.01:xp2;
x3=xp2:0.01:L;

A1=C*(xp1^2-x1.^2);
A2=C*(-x2.^2+(xp1+xp2)*x2-xp1*xp2);
A3=C*((L-xp2)^2-(L-x3).^2);

plot(x1,A1,x2,A2,x3,A3);
sum([A1,A2,A3])


figure(3)
%%with global unbinding/degredation

xp=0.3; %plasmid position
xl=0:0.01:xp;
xr=xp:0.01:L;

kp=1;%production rate
d=1;%global degradation rate
lengthscale=0.2;
d=2*D/(lengthscale^2);
%sqrt(2*D/d)

Al=kp/d-kp/d*(exp(sqrt(d/D)*xl)+exp(-sqrt(d/D)*xl))/(exp(sqrt(d/D)*xp)+exp(-sqrt(d/D)*xp));
Ar=kp/d-kp/d*(exp(sqrt(d/D)*(xr-L))+exp(-sqrt(d/D)*(xr-L)))/(exp(sqrt(d/D)*(xp-L))+exp(-sqrt(d/D)*(xp-L)));
m=max([Al,Ar]);
Al=Al/m;
Ar=Ar/m;

subplot(3,1,1)

plot(xl,Al,xr,Ar)
ylim([0,1.1])
ylabel('Concentration/mass')

% subplot(3,1,2)
figure(4)
dflux=[];
Atot=[];
for xp=0:0.01:1
xl=0:0.01:xp;
xr=xp:0.01:L;
Al=kp/d-kp/d*(exp(sqrt(d/D)*xl)+exp(-sqrt(d/D)*xl))/(exp(sqrt(d/D)*xp)+exp(-sqrt(d/D)*xp));
Ar=kp/d-kp/d*(exp(sqrt(d/D)*(xr-L))+exp(-sqrt(d/D)*(xr-L)))/(exp(sqrt(d/D)*(xp-L))+exp(-sqrt(d/D)*(xp-L)));
fluxL=D*abs(gradient(Al))/0.01;
fluxR=D*abs(gradient(Ar))/0.01;
dflux=[dflux fluxR(1)-fluxL(end)];
Atot=[Atot sum([Al,Ar])];

end
x = 0:0.01:1;
x = x-0.5;
plot(x,dflux,'LineWidth',2)
% ylabel('Flux differential')
% xlabel('x')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% subplot(3,1,3)
% plot(0:0.01:1,Atot)
% ylabel('Total mass')