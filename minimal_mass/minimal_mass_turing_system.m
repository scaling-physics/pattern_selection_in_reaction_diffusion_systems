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

k = 0.4327;

%changing k
% k = 0.1:0.01:10;
b = k.^2/L^2*D/gamma;
% b = 0.0039;

delta = b*gamma;
% k = L*sqrt(delta/D);
beta = a*gamma/c^2;

%% Plasmid parameters mapping
n = 6;
% changing n
n=1:200;
gp = coth(k./(2.*n));
rho = 6.*sqrt(Dv)/L.*(delta).*sqrt(delta + gamma)/beta;
muprime1 = c.*delta./(gp.*k);
muprimep = muprime1 + sqrt(muprime1.^2 - (2.*delta)./(gp.*k).*rho);
muprimem = muprime1 - sqrt(muprime1.^2 - (2.*delta)./(gp.*k).*rho);
% 
% M1 = c - n./(delta).*muprimep;
% M2 = c - n./(delta).*muprimem;
% 
% % A = 1./M1(~imag(M1));
% % find(A == max(A))
% %% Plot
% figure(1)
% hold on
% plot(n(~imag(M1)),1./M1(~imag(M1))/c)
% plot(k(~imag(M1)),M1(~imag(M1))/c)
% plot(n,M1(~imag(M1)))
% hold on
% plot(n(~imag(M2)),1./M2(~imag(M2)))
% end
% xlabel('n ','FontSize', 18);
% xlabel('\kappa ','FontSize', 18);
% ylabel('$\frac{1}{Lc}\int u dx$','FontSize', 18,'interpreter','latex');

% legends
% k1 = 0.4327;
% n1= 1;
% legend(['n=' num2str(n1)], ['n=' num2str(n1+1)], ['n=' num2str(n1+2)], ['n=' num2str(n1+3)], ['n=' num2str(n1+4)], ['n=' num2str(n1+5)],['n=' num2str(n1+6)]);
% legend(['k=' num2str(k1)], ['k=' num2str(k1*2)], ['k=' num2str(k1*4)], ['k=' num2str(k1*8)], ['k=' num2str(k1*12)], ['k=' num2str(k1*16)],['k=' num2str(k1*20)]);
% legend(l);
% end
