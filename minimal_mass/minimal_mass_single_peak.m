%% SMC parameters
L = 2;
dx = 0.01;
% beta = 1.5e-4;/
gamma = 3.6;
D = 0.3;
s = 1;
Dv = 0.012*s;

c = 300;
a = 3.75;

% k = 0.4327;

%changing k
% b = k.^2/L^2*D/gamma;
b = 0.0039;
delta = b*gamma;
k = L*sqrt(delta/D);
beta = a*gamma/c^2;
rho = 6.*sqrt(Dv)/L.*(delta).*sqrt(delta + gamma)/beta;

%% Minimal mass single peak
x1 = -L/2:dx:L/2;
x = -L/2:dx:L/2;
gprime = k/2.*(cosh(2*k.*x1/L)/sinh(k)+coth(k));
bdash = c*delta./gprime;
rhoprimep = 0.5*bdash.*(1 + sqrt(1 - 4*rho.*gprime./(delta*c^2)));



%% Utotal
load('pattern.mat');
load('utotal.mat');

figure(1)
gprime0 = k/2.*(cosh(k*x/L)+cosh(k*(abs(x)-L)/L))/sinh(k);
% plot(x/L,gprime0)
gprime00 = k/2*coth(k/2);
rhoprime0 = 0.5*c*delta./gprime00*(1 + sqrt(1 - 4*rho*gprime00/(delta*c^2)));
up =   c - rhoprime0/delta*gprime0;

plot(x/L,up/c)
hold on

plot(x/L,u1(end,:)/(c))
xlabel('x/L','FontSize', 18);
ylabel('u(x)/c','FontSize', 18);
%%  u(x)
figure(2)
plot(x/L,u1(end,:)/c)

I = get_last_monotonic(utotal,7)+20;
utotal_a = c - rhoprimep/delta;
 
figure(3)
plot(x1/L, utotal_a/(c))
hold on
plot(centroid(I:end)/L, utotal(I:end),'LineWidth',2,'color','r');
plot(-flip(centroid(I:end))/L, flip(utotal(I:end)),'LineWidth',2,'color','r');

xlabel('x/L');
ylabel('$\frac{\int u}{Lc}$','interpreter','latex')

%% Flux
deltaJ = -L/2*rhoprimep.*sinh(2*k*x1/L)/sinh(k);
I=1;
vtotal_a = c-utotal_a;
theta = c*D/L;
load('flux.mat')
figure(4)
plot(x1/L,deltaJ)
set(gca,'YAxisLocation','Origin','XAxisLocation','Origin')
hold on
plot(centroid(I:end)/L, flux(I:end,1),'Color','b');
plot(-centroid(I:end)/L,-flux(I:end,1),'Color','b');

% xlabel('x1/L','FontSize', 18);
% ylabel('\Delta J','FontSize', 18);

