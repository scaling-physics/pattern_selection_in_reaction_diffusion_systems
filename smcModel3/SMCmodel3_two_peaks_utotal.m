
function [peak1,peak2,utotal] = SMCmodel3_two_peaks_utotal(shift,shift1)
% function SMCmodel3_two_peaks_utotal
% addpath('/Users/srik/Documents/MATLAB/export_fig-master');
%%%%%   in units of nM
%%%% u(u+v)^2
s = 1; %scaling up
L = 4;
m = 0;
x = 0:0.01:L;
t = 0:100:60*60*5;

% k1 = s*1.5e-4;
k_minus1 = s*3.6;
Du = s*0.3;
Dv = s*0.012*0.1;
% koff = 0.001*log(2)/50;


C=300;



%old definitions
% b=(koff+k_minus1)/k_minus1;
% v0=C*a/(a+b);
% u0=C-v0;
%new defintions
b = 0.0039;
% k = 0.0139;
k = L*sqrt(b*k_minus1/Dv);
% b = k^2/(L^2)*Dv/(k_minus1);
a = 3.75;
koff = k_minus1*b;
k1 = a*k_minus1/(C^2);
kon = C*koff;
% a = k1/k_minus1*C^2;
% b = koff/k_minus1;

Gamma = k_minus1*L^2/Dv;

v0 = C*a/(a+b+1);
u0 = C-v0;


data = load('pattern_two_peaks_spike.mat');
% shift = 140;
% shift1 = 120;

u20 = data.u2(end,1:round(end/2)); 

u21 = u20(end,[end-shift:end,1:end-shift-1]);
u20 = data.u2(end,round(end/2):end); 
u22 = u20(end,[end-shift1:end,1:end-shift1-1]);
data.u2 = [u21 u22(end,1:end-1)];
% v00 = v0*exp(-25*(x).^2);
% v1 = v00([end-shift:end,1:end-shift-1]);
% v2 = v00([end-shift1:end,1:end-shift1-1]);
% v0 = v1+v2;
% data.u2 = data.u2(end,[end-shift:end,1:end-shift-1]);
% data.u1 = data.u1(end,[end-shift:end,1:end-shift-1]);
v0 = data.u2;
u0 = data.u1(end,:);
IC = [u0;v0];

% 
% rr1=randn(size(x));
% rr2=randn(size(x));
% rr2=sum(rr1)*rr1/sum(rr2);
% IC = [u0*(1+0.1*rr1);v0-0.1*u0*rr2];

options=odeset('RelTol',1e-6,'AbsTol',1e-12);
sol = pdepe(m,@pdes,@ic,@bc,x,t,options);
u1 = sol(:,:,1);
u2 = sol(:,:,2);


% save('pattern_two_peaks_spike.mat','u1','u2')

%%
% figure(1)

peak1 = [];
peak2 = [];

for i=1:length(t)
[Peaklist, PFWHH] = mspeaks(x,u2(i,:),'HeightFilter',250,'PeakLocation', 0.4);%position of peaks
if length(Peaklist(:,1))==2
peak1 = [peak1 Peaklist(1,1)];
peak2 = [peak2 Peaklist(2,1)];
end
end
peak1 = peak1'-L/2;
peak2 = peak2'-L/2;
utotal = sum(u1,2)/(L/0.01);
if utotal
% I = get_last_monotonic(gradient(utotal),7);
I=1;
peak1= peak1(I:end)/L;
peak2 = peak2(I:end)/L;
utotal = utotal(I:end)/C;
end
scatter(peak1,peak2)
set(gca,'XDir','reverse')
% hold on
% scatter(peak1(I:end)/L,peak2(I:end)/L,10,utotal(I:length(peak1))/C,'filled');
% scatter(-peak1(I:end)/L,-peak2(I:end)/L,10,utotal(I:length(peak1))/C,'filled');
% 
% colorbar()
xlim([-0.5 0])
ylim([0 0.5])
% 
% m = [peak1(I:end)/L; peak2(I:end)/L; utotal(I:length(peak1))/C];




%%
% figure(3)
% %[~,I]=max(u2,[],2);
% newdata=load('flux_balance.mat','centroid','flux');
% 
% centroid=diag(1./sum(u2,2))*u2*x'; %position of centroid
% 
% % Plot fit with data.
% centroid = centroid-(L/2);
% plot(t/60,centroid);
% plot(t/60,newdata.centroid);
% 
% I=get_last_monotonic(gradient(centroid),7);
% % I1=get_last_monotonic(gradient(newdata.centroid),7);
% %% Flux across the peak
% % flux = diag(1./sum(u2,2))*trapz(gradient(u1).*u2)
% num = gradient(u1,2).*u2;
% num = trapz(num,2);
% flux = num./trapz(u2,2);
% 
% scatter(centroid(I:end), flux(I:end),10,'filled','blue');
% 
% hold on
% scatter(newdata.centroid(I:end), newdata.flux(I:end),10,'filled','blue');
% f=robustfit(centroid(I:end),flux(I:end));
% k=robustfit(newdata.centroid(I1:end),newdata.flux(I1:end));
% plot(centroid(I:end),f(1)+f(2)*centroid(I:end),'r','LineWidth',1.5);
% plot(newdata.centroid(I1:end),k(1)+k(2)*newdata.centroid(I:end),'r','LineWidth',1.5);
% hold off
% % xlim([-0.5,0.5])
% % ylim([-max(flux), max(flux)])
% ax=gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% xlabel('Position');
% ylabel('Flux');
%%

% U=u2(end,:);
% U=[0,U,0];
% xx=[-0.05,x,L+0.05];
% Y=mspeaks(xx,U,'heightfilter',300);



%%
% 
% figure(2);
% clf;
% % subplot(2,1,1);
% % imagesc(t/60,x-L/2,u1');
% % %set(p1,'LineStyle','none');
% % title('u(x,t)')
% % ylabel('Position (micro m)')
% % xlabel('Time (min)')
% % zlim([0,5]);
% 
% % subplot(2,1,2);
% imagesc(t/60,(x-L/2)./L,u2');
% % set(p2,'LineStyle','none')
% title('v(x,t)')
% ylabel('x/L')
% xlabel('Time (min)')
% zlim([0,15]);
% export_fig('/Users/srik/Documents/Pattern_selection_images/Figure 2/peak_kymograph.pdf');
% 
% %%
% 
% 
% figure(3);
% clf;
% % subplot(4,1,1);
% plot(x,u1(end,:),x,u2(end,:))
% title('Solution at t = end');
% xlabel('Distance x');
% % ylim([0,2*C]);
% 
% figure(3)
% % subplot(4,1,2)
% for n=0:40
% an(n+1)=2*trapz(x,u2(end,:).*cos(pi*n*x/L)/1000)/L;
% end
% plot(0:15,an(1:16),'-x','LineWidth',2)
% line([0,15],[0,0])
% xlabel('mode')
% ylabel('amplitude')
% a = an(2:16);
% [val,idx] = max(abs(a))
% 
% figure(4)
% % subplot(4,1,3);
% for i=1:length(t)
% s1(i)=trapz(x,u1(i,:))/L;
% % s2(i)=trapz(x,u2(i,:))/L;
% end
% plot(t/60,s1)
% % plot(t/60,s1,t/60,s2,t/60,s1+s2);
% % title('U');
% ylabel('$\bar{u}/L$ ','interpreter','latex');
% xlabel('Time (min)');
% %ylim([0,2*C]);
% 
% % subplot(4,1,4);
% % [t1,Y]=ode23s(@odes,t,[u0,v0],options);
% % plot(t1/60,Y(:,1),t1/60,Y(:,2),t1/60,Y(:,1)+Y(:,2))
% % title('Total concentrations (0-d model)');
% % ylabel('nM');
% % xlabel('Time');
% % 
% % ylim([0,2*C]);
% 
% 
% 
% 
% Centroid rate

%[~,I]=max(u2,[],2);
% centroid=diag(1./sum(u2,2))*u2*x'; %position of centroid
% 
% 
% % Plot fit with data.
% figure(5)
% centroid = centroid-(L/2);
% % plot(t/60,centroid,'LineWidth',2);
% 
% I = get_last_monotonic(gradient(centroid),7);
% 
% % Flux across the peak
% % flux = diag(1./sum(u2,2))*trapz(gradient(u1).*u2)
% num = gradient(u1,2).*u2;
% num = trapz(x,num,2);
% flux = Du*num./trapz(x,u2,2);
% % I = get_last_monotonic(gradient(centroid),7);
% dx = 0.01;
% theta = C*Du/L;
% 
% 
% plot(centroid(I:end)/L, flux(I:end,1),'LineWidth',2,'Color','b');
% hold on
% plot(-centroid(I:end)/L,-flux(I:end,1),'Color','b','LineWidth',2);
% set(gca,'YAxisLocation','Origin','XAxisLocation','Origin')
% xlabel('x1/L');
% ylabel('$\theta^{-1}\Delta J$','interpreter','latex');
% xlim([-1/2 1/2])
% save('flux.mat','flux','centroid')
% lgd = legend('Simulation','Fit', 'Location', 'north east');

% figure(6)
% % Centroid rate fit
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
% plot(t(I:end)/60, feval(fitresult,  t(I:end)'-t(I)),'LineWidth',2)
% hold off;
% % ylim([-0.3 0.3]);
% xlabel('time (min)');
% ylabel('x_1/L')
% 
% 
% 
% figure(7)
% I = get_last_monotonic(gradient(centroid),7);
% % plot(feval(fitresult,  t(I:end)'-t(I))/L,gradient(feval(fitresult,  t(I:end)'-t(I))),'Color','b');
% % hold on
% % plot(-feval(fitresult,  t(I:end)'-t(I))/L,-(gradient(feval(fitresult,  t(I:end)'-t(I)))),'Color','b');
% plot(centroid(I:end)/L,gradient(centroid(I:end)),'Color','b');
% hold on
% plot(-centroid(I:end)/L,-(gradient(centroid(I:end))),'Color','b');
% set(gca,'YAxisLocation','Origin','XAxisLocation','Origin')

% xlabel('x_1/L');
% ylabel('Velocity')

%% mode 4 fit
% figure(3)
% peak1 = [];
% peak2 = [];
% 
% for i=1:length(t)
% [Peaklist, PFWHH] = mspeaks(x,u2(i,:),'HeightFilter',250,'PeakLocation', 0.4);%position of peaks
% if length(Peaklist(:,1))==2
% peak1 = [peak1 Peaklist(1,1)];
% peak2 = [peak2 Peaklist(2,1)];
% end
% end
%  
% 
% peak1 = peak1';
% peak2 = peak2';
% I1 = length(t)-length(peak1)+1;
% d1 = peak1 + peak2 - L;
% d2 = peak1 - peak2 + L/2;
% w = 1./t;
% w = w(I1:end);
% ft = fittype('a*exp(-b*t)', 'independent', 't', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Weight = w;
% opts.Lower = [0 0];
% opts.StartPoint = [koff d1(1)];
% opts.Upper = [Inf L];
% % % 
% % % Fit model to data.
% 
% [fitresult1, gof1] = fit(t(I1:end)'-t(I1), d1, ft, opts);
% rate1=coeffvalues(fitresult1);
% % rate1 = rate1(2)
% rate1
% I2 = length(t)-length(peak2)+1;
% w = 1./t;
% w = w(I2:end);
% 
% ft = fittype('a*exp(-b*t)', 'independent', 't', 'dependent', 'y' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Weight = w;
% opts.Lower = [0 0];
% opts.StartPoint = [koff d2(1)];
% opts.Upper = [Inf L];
% % 
% [fitresult2, gof2] = fit(t(I2:end)'-t(I2),d2, ft, opts);
% rate2=coeffvalues(fitresult2);
% rate2
% % rate2 = rate2(2)
% subplot(2,1,1);
% plot(t(I1:end)/60,d1);
% hold on
% plot(t(I1:end)/60, feval(fitresult1,  t(I1:end)'-t(I1)))
% hold off
% ylabel('d_1 + d_2')
% subplot(2,1,2);
% plot(t(I2:end)/60,d2);
% hold on
% plot(t(I2:end)/60, feval(fitresult2,  t(I2:end)'-t(I2)))
% hold off;
% xlabel('time (min)');
% ylabel('d_1 - d_2')

%%
%phase portrait

% y1max=400;
% y2max=400;
% %y1=0:1e-8/unit:y1max;
% %y2=0:1e-8/unit:y2max;
% y1 = linspace(0,y1max,50);
% y2 = linspace(0,y2max,50);
% 
% [x,y] = meshgrid(y1,y2);
% u = zeros(size(x));
% v = zeros(size(x));
% 
% for i = 1:numel(x)
%     [~,~,Yprime] = pdes(0,0,[x(i); y(i)],0);
%     u(i) = Yprime(1);
%     v(i) = Yprime(2);
% end
% 
% for i = 1:numel(x)
% Vmod = sqrt(u(i)^2 + v(i)^2);
% u(i) = u(i)/Vmod;
% v(i) = v(i)/Vmod;
% end
% 
% 
% figure(6);
% clf;
% quiver(y,x,v,u,'r');
% hold on;
% %axis tight equal;
% 
% null1=@(v,u) kon-koff*u-koff*v;
% null2=@(v,u) k1*u.*(u+v).^2-k_minus1*v-koff*v;
% p6=ezplot(null1,[0,y2max,0,y1max]);
% p7=ezplot(null2,[0,y2max,0,y1max]);
% set(p6,'Color','k','Linewidth',2);
% set(p7,'Color','k','Linewidth',2);
% xlabel('v (nM)')
% ylabel('u (nM)')
% 
% plot(Y(:,2),Y(:,1),'Linewidth',3);
% plot(Y(1,2),Y(1,1),'bo') % starting point
% plot(Y(end,2),Y(end,1),'ks') % ending point
% 
% plot(s2,s1,'Linewidth',3,'Color','g');
% plot(s2(1),s1(1),'bo') % starting point
% plot(s2(end),s1(end),'ks') % ending point
% 
% hold off;
% xlim([0,400])


%%
%Dispersion Relation

% %Jacobian
% fu=-k1*(u0+v0)^2-2*k1*u0*(u0+v0)-koff;
% fv=-2*k1*u0*(u0+v0)+k_minus1;
% gu=k1*(u0+v0)^2+2*k1*u0*(u0+v0);
% gv=2*k1*u0*(u0+v0)-k_minus1-koff;
% 
% 
% function [RE, IM]=A_e_values(k2)
%      AA=[fu-Du*k2, fv;
%      gu, gv-Dv*k2];
%     lambda=eig(AA);
%     [RE,I]=max(real(lambda));
%     IM=imag(eig(I));
% end
% 
% nrange=[0:1:4*L];
% krange=0:0.05:pi*nrange(end)/L;
% 
% for i=1:length(nrange)
%     %k=pi*n/L;
%     k=pi*nrange(i)/L;
%     [A(i),A_im(i)]=A_e_values(k.^2);
% end
% 
% 
% 
% %%alternative: explicit formula for the greater eigenvalue
% %Eig=@(k2) ((fu+gv-k2*(Du+Dv))+sqrt((fu+gv-k2*(Du+Dv)).^2-4*(fu*gv-gu*fv-k2*(fu*Dv+gv*Du)+k2.^2*Du*Dv)   )  )/2;
% 
% figure(8);
% plot(nrange,A,'-o',nrange,A_im, 'LineWidth',2)
% % title('Dispersion Relation','FontSize',12);
% xlabel('n')
% ylabel('Growth rate')
% set(gca,'FontName', 'Helvetica','FontSize', 14)
% set(gca,'color','none');
% % export_fig('/Users/srik/Documents/Pattern_selection_images/Figure 1/linear_dispersion_relation.pdf','-transparent');
% % a = get(gca,'XTickLabel');
% % 
% % set(gca,'XTickLabel',a,'fontsize',25)
% % % lgd = legend('Centroid Position','Exponential Fit', 'Location', 'north east');
% % set(gca,'Color', 'none');
% 
% % savefig('/Users/subraman/Documents/MATLAB/SMC/data/poster_images/disp_relation.eps');
% 
% 
% %%b
% %Turing space
% 
% figure(9)
% % cond=@(d,b,a) d*b.*(a-b)./(a+b)+(-(a.^2+b.^2+4*a.*b)./(a+b)+1)-2*sqrt(d*(b-1).*(a+b)); %original definition of b
% cond=@(d,b,a) d*(b+1).*(a-b-1)./(a+b+1)+(-(a.^2+b.^2+4*a.*b+3.*a+b)./(a+b+1))-2*sqrt(d*(b).*(a+b+1)); %new b definition
% 
% d=Du/Dv;
% h = ezplot(@(b,a) cond(d,b,a),[0,5,0,30]);
% set(h,'LineWidth',2);
% 
% hold on
% % plot(b,a, 'LineWidth', 1);
% hold off;
% xlim([0,0.25])
% ylim([0,25])

%% Utotal
% figure(8)
% 
% utotal = sum(u1,2)/(L/0.01)/C; %position of centroid
% % utotal = trapz(x,u1,2)/(L*C);
% % Plot fit with data.
% % centroid = centroid-(L/2);
% 
% centroid = diag(1./sum(u2,2))*u2*x'-L/2;
% I = get_last_monotonic(utotal,7);
% plot(centroid(I:end)/L,utotal(I:end),'Color',[0    0.4470    0.7410]);
% hold on
% plot(-centroid(I:end)/L,flip(utotal(I:end),2),'Color',[0    0.4470    0.7410]);
% xlabel('x_1/L')
% ylabel('$\frac{1}{Lc}\int u$', 'interpreter', 'latex')
% save('utotal_spike.mat','utotal','centroid')

%% Utotal two peaks
% figure(9)
% utotal = sum(u1,2)/(L/0.01);
% peak1 = [];
% peak2 = [];
% 
% for i=1:length(t)
% [Peaklist, PFWHH] = mspeaks(x,u2(i,:),'HeightFilter',250,'PeakLocation', 0.4);%position of peaks
% if length(Peaklist(:,1))==2
% peak1 = [peak1 Peaklist(1,1)];
% peak2 = [peak2 Peaklist(2,1)];
% end
% end
% peak1 = peak1'-L/2;
% peak2 = peak2'-L/2;
% I = get_last_monotonic(gradient(utotal),7);
% I = I+15;
% 
% if(abs(abs(peak1(end)-peak2(end)) - L/2) < 0.1)
% scatter3(peak1(I:end)/L,peak2(I:end)/L,utotal(I:length(peak1))/C,10,'filled')
% end
% hold on
% xlabel('x_1')
% ylabel('x_2')
% zlabel('$\frac{1}{Lc}\int u$', 'interpreter', 'latex')
% xlim([-1/2 0])
% ylim([0 1/2])
%%
% figure(5)
% eta = (Du/Dv)*u1(end,:) + u2(end,:);
% plot(x,eta)
% 
% 
% figure(6)
% function dydt = f(x,y)
% 
% dydt = [y(2); ((b/Du)*y(1)-b)];
% end
% [x,y] = ode45(@f,[0 1],[1620 ; 0]);
% 
% plot(x,y(:,1),'-o')
%%


function out=odes(~,u)
    [~,~,out]=pdes(0,0,u,0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------
function [c,f,s] = pdes(x,t,u,DuDx)

c = [1; 1];
f = [Du;Dv] .* DuDx; 
F = k1*u(1)*(u(1)+u(2))^2-k_minus1*u(2);%
s = [-F+kon-koff*u(1);F-koff*u(2)]; 
end
% --------------------------------------------------------------
function out = ic(xx)
out = IC(:,x==xx);
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bc(xl,ul,xr,ur,t)
pl = [0; 0]; %negative=flux out, positive= flux in
ql = [1; 1]; 
pr = [0; 0]; % negative=flux in, positive=flux out
qr = [1; 1];
end
% function [pl,ql,pr,qr] = bc(xl,ul,xr,ur,t)
% pl = [ul(1)-ur(1); ul(2)-ur(2)]; %negative=flux out, positive= flux in
% ql = [0; 0]; 
% pr = [0; 0];
% qr = [1; 1];
% end




end