function SMCmodel3_spike

filename = '/Users/srik/Documents/MATLAB/flow/data/raw_data/spike_solution.mat';
set(0,'defaultAxesFontSize',20)
%Parameters
m = 0;
L = 2;
x = -L/2:0.01:L/2;
% x = x./L;

t = 0:60:60*60*10;
C = 300;
k1 = 1.5e-4;
k_minus1 = 3.6;
koff = log(2)/50;
% koff = 0.;
a = k1/k_minus1*C^2;
b = koff/k_minus1;
v0 = a/(a+b+1);
u0 = 1-v0;

Du = 0.3;
Dv = 0.012;
D = Du./(k_minus1*L^2);
% D = 0.08;
% epsilon = sqrt(Dv./(k_minus1*L^2));
epsilon = 0.01;

A = a/epsilon;
k = sqrt(b/D);



%Initial condition with peak
data = load(filename);

shift = 40;

% data.u2 = data.u2(end,[end-shift:end,1:end-shift-1]);
% data.u1 = data.u1(end,[end-shift:end,1:end-shift-1]);
% data.u2 = data.u2(end,[end-shift:end,1:end-shift-1]);
% data.u1 = data.u1(end,[end-shift:end,1:end-shift-1]);

% v0 = data.u2;
% u0 = data.u1;
% IC = [u0;v0];



% Initial conditions
rr1=randn(size(x));
rr2=randn(size(x));
rr2=sum(rr1)*rr1/sum(rr2);
IC = [u0*(1+0.1*rr1);v0-0.1*u0*rr2];

options=odeset('RelTol',1e-6,'AbsTol',1e-12);
sol = pdepe(0,@pdes,@ic,@bc,x,t,options);
u1 = sol(:,:,1);
u2 = sol(:,:,2);

% l=1;
% %save plots
% for k=1:100:length(t)
%     clf
%     figure(4);
%     plot(x,u2(k,:),'Color',[255, 89, 105]./256,'LineWidth',4);
% %     hold on
% %     plot(x,u1(k,:),'Color',[82, 203, 190]./256,'LineWidth',4);
% %   ylabel('v,u');
% %     xlabel('x');
%     ylim([0,20]);
%     axis off ; 
%     set(gca,'yticklabel',[],'xticklabel',[],'XTick',[],'YTick',[],'box','off');
%     set(gca,'Color', 'none');
%     export_fig(sprintf('/Users/subraman/Documents/MATLAB/SMC/data/animation_figures/2/figure%03d.png', l),'-transparent');
%     l=round(k/100);
% end
%% Data comparisons
ubar_final = 6*sqrt(b+1)/A*k.*cosh(k/2)./sinh(k).*cosh(k/2);
v_data = trapz(x,u2(end,:))
v_analytical = 6*sqrt(b+1)/(A*ubar_final)
rate_a = epsilon^2*6*b*sqrt(b+1)/(D*A*ubar_final^2);
n = C/koff;
% %% Plot
figure(1);
clf;
% subplot(4,1,1);
plot(x,u2(end,:),x,u1(end,:))
% title('Solution at t = end');
xlabel('Distance x');
ylabel('u,v')
export_fig('/Users/srik/Documents/Presentation_images/spike_solution.pdf','-transparent')
%ylim([0,2*C]);

% figure(2) 
% clf;
% subplot(2,1,1);
% imagesc(t/60,x,u1');
% %set(p1,'LineStyle','none');
% title('u(x,t)')
% ylabel('Position (micro m)')
% xlabel('Time (min)')
% %zlim([0,5]);

% subplot(2,1,2);
figure(3)
imagesc(t,x,u2');
% hold on
% imagesc(t,x,u1');
%set(p2,'LineStyle','none')
title('v(x,t)')
ylabel('Position (micro m)')
xlabel('Time (s)')
%zlim([0,15]);

%%save pattern
% save(filename, 'u1','u2')

%% Subspace constant
ubar = trapz(x,u1(end,:));
vbar = trapz(x,u2(end,:));
C1 = D*ubar + epsilon^2*vbar
C1_analytical = 3*sqrt(D/2*a)*epsilon 

%% Numerical comparison of approximations
figure(5)
eta = D*u1(end,:) + epsilon^2*u2(end,:);
vout = u2(end,1);
plot(x,eta/D)
hold on
plot(x,u1(end,:))

% function dydt = f(x,y)
% 
% dydt = [y(2); ((b/D)*y(1)-b)];
% end
% [x1,y] = ode45(@f,[-0.5 0],[eta(1) ; 0]);
% 
% hold on
% plot(x1,y(:,1),'-o','Color','r')
% [x1,z] = ode45(@f,[0.5 0],[eta(end) ; 0]);
% plot(x1,z(:,1),'-o','Color','r')
% 
% xlabel('Position');
% ylabel('\eta');
ylim([0 0.1]);



% % figure(6)
plot(x,u1(end,:)); 
function dydt = f1(x,y)
dydt = [y(2); b/D*(y(1)-1)];
end
[x2,y] = ode45(@f1,[-0.5 0],[u1(end,1) ; 0]);

hold on
plot(x2,y(:,1),'-o','Color','r')
[x2,z] = ode45(@f1,[0.5 0],[u1(end,end) ; 0]);
plot(x2,z(:,1),'-o','Color','r')
xlabel('Position');
ylabel('\eta,u');
% 
% 
% figure(7)
% plot(x,u2(end,:));
% figure(8)
% y = 3*(b+1)/(2*A*ubar)*sech(sqrt(b+1)/2.*(x-0.5)/epsilon);
% plot(x,y);
%% Series expansion estimate
eta1 = trapz(x,eta);
v0 = u2(end,51);
c1 = ((a/D)*epsilon^2*v0^2+(b+1)-(a/D)*eta1*v0)/(2*epsilon^2);
c2 = (3*(a/D)*epsilon^2*v0^2+(b+1)-2*(a/D)*eta1*v0)/(12*epsilon^2);
v_pade_series = v0.*(1+(c1-c2).*x.^2)./(1-c2.*x.^2);
% v_series = v0.*(1+c1*x.^2+c1*c2*x.^4);
figure(10)
plot(x,v_pade_series)
% hold on
% plot(x,v_series)

%% Flux across the peak
% figure
% %[~,I]=max(u2,[],2);
% centroid=diag(1./sum(u2,2))*u2*x'; %position of centroid
% 
% % Plot fit with data.
% % centroid = centroid-(L/2);
% % plot(t/60,centroid);
% 
% I=get_last_monotonic(gradient(centroid),7);
% 
% 
% [uxgrad,~] = gradient(u1,0.001);
% num = uxgrad.*u2;
% num = trapz(x,num,2);
% flux = num./trapz(x,u2,2);
% % flux = num; %definition without normalization
% centroid = centroid(I:end);
% flux = flux(I:end);
% 
% 
% scatter(centroid, flux,10,'filled');
% f = robustfit(centroid,flux);
% 
% x1 = 0:0.01:0.5;
% k = sqrt(b/D);
% m = 6*sqrt(b+1)/A*k.*cosh(k.*x1)./sinh(k).*cosh(k.*(1-x1));
% c1 = -0.5*k*sinh(k.*(2.*x1-1))./(cosh(k.*x1).*cosh(k.*(1-x1)));
% 
% 
% hold on
% plot(x1,c1);
% % plot(centroid,f(1)+f(2)*centroid,'r','LineWidth',1.5);
% lgd = legend('Simulation','Fit', 'Location', 'north east');
% xlabel('Centroid');
% ylabel('Flux');

%% save file
% file = '/Users/subraman/Documents/MATLAB/SMC/data/raw_data/flux_on_spike_1.mat';
% % save(file)
% m = matfile(file,'Writable',true);
% % % 
% m.centroid = [m.centroid; centroid];
% m.flux = [m.flux; flux]


% -------------------------------------------------------------
% Centroid Rate
% figure(4)
% %[~,I]=max(u2,[],2);
% centroid=diag(1./sum(u2,2))*u2*x'; %position of centroid
% 
% % Plot fit with data.
% centroid = centroid-(L/2);
% plot(t/60,centroid);
% 
% I=get_last_monotonic(gradient(centroid),7);
% %% Flux across the peak
% % flux = diag(1./sum(u2,2))*trapz(gradient(u1).*u2)
% num = gradient(u1,2).*u2;
% num = trapz(num,2);
% flux = num./trapz(u2,2);
% scatter(centroid(I:end), flux(I:end),10,'filled');
% f=robustfit(centroid(I:end),flux(I:end));
%  
% hold on
% plot(centroid(I:end),f(1)+f(2)*centroid(I:end),'r','LineWidth',1.5);
% xlabel('Centroid');
% ylabel('Flux');
% lgd = legend('Simulation','Fit', 'Location', 'north east');

%% Slowly moving pulse

rate_a = epsilon^2*b*A/(3*D);
m = epsilon^2*A/(6*D*k_minus1);


%% Centroid rate fit
% figure(3)
% L=1;
% centroid=diag(1./sum(u2,2))*u2*x'; %position of centroid
% 
% % Plot fit with data.
% centroid = centroid-(L/2);
% plot(t/60,centroid);
% 
% I=get_last_monotonic(gradient(centroid),7);
% 
% ft = fittype('a*exp(-b*t)', 'independent', 't', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Lower = [-L/2 0];
% opts.StartPoint = [1 koff];
% opts.Upper = [L/2 1];
% 
% % Fit model to data.
% [fitresult, gof] = fit( t(I:end)'-t(I), centroid(I:end), ft, opts );
% rate=coeffvalues(fitresult);
% rate
% hold on;
% 
% plot(t(I:end)/60, feval(fitresult,  t(I:end)'-t(I)));
% hold off;
% xlabel('Time (min)');
% ylabel('Peak Position');
% lgd = legend('Simulation','Fit', 'Location', 'north east');

%% Velocity of Spikes
% figure
% 
% position_sim = centroid(I:end);
% velocity_sim =  gradient(feval(fitresult,  t(I:end)'-t(I)));         
% plot(position_sim,velocity_sim);
%% Saving data
% file = '/Users/subraman/Documents/MATLAB/SMC/data/raw_data/spike_velocity_analysis_low_k.mat';
% % save(file)
% % 
% m = matfile(file,'Writable',true);                                                                                                                                                                                    
%  
% m.position_sim = [m.position_sim;position_sim];
% m.velocity_sim = [m.velocity_sim;velocity_sim];


%% --------------------------------------------------------------
function [c,f,s] = pdes(~,~,u,DuDx)

c = [1; 1];
f = [D;epsilon^2] .* DuDx; 
F = a*u(1)*(u(1)+u(2))^2-u(2);
s = [-F+b-b*u(1);F-b*u(2)]; 
end
% --------------------------------------------------------------
function out = ic(xx)

out = IC(:,x==xx);
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bc(xl,ul,xr,ur,t)
% pl = [ul(1)-1;ul(2)]; 
% ql = [0; 0]; 
% pr = [ur(1)-1;ur(2)];
% qr = [0; 0];
pl = [0; 0]; %negative=flux out, positive= flux in
ql = [1; 1]; 
pr = [0; 0]; % negative=flux in, positive=flux out
qr = [1; 1];
end


end