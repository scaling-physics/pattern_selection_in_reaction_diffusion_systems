function schnakenberg_model

m = 0;
L= 4;
x = 0:0.01:L;
t = 0:1:600;

C1 = 1;
C2 = 1;

Da = 0.01;
Db = 1;


k2 = 0.2;
k3 = 2;
k1 = C1*k2;
k4 = C1*k2;

u0 = C1;
v0 = k2/(k3*C1);
% u0 = a+b;
% v0 = b/(a+b)^2;


a1 = k3/k2
gamma = k2*L^2/Da
% b = (k4/k2)*(k3/k2)^0.5

rr1=randn(size(x));
rr2=randn(size(x));
rr2=sum(rr1)*rr1/sum(rr2);
IC = [u0*(1+0.1*rr1);v0*(1-0.1*rr2)];



options=odeset('RelTol',1e-6,'AbsTol',1e-12);
sol = pdepe(0,@pdes,@ic,@bc,x,t,options);
u1 = sol(:,:,1);
u2 = sol(:,:,2);



%% Plot
figure(1);
clf;
subplot(4,1,1);
plot(x,u2(end,:),x,u1(end,:))
title('Solution at t = end');
xlabel('Distance x');
%ylim([0,2*C]);

figure(2); 
clf;
subplot(2,1,1);
imagesc(t/60,x,u2');
%set(p1,'LineStyle','none');
title('u(x,t)')
ylabel('Position (micro m)')
xlabel('Time (min)')
%zlim([0,5]);

subplot(2,1,2);
imagesc(t/60,x,u1');
%set(p2,'LineStyle','none')
title('v(x,t)')
ylabel('Position (micro m)')
xlabel('Time (min)')
%zlim([0,15]);



%% Centroid rate fit
figure(3)

centroid=diag(1./sum(u2,2))*u2*x'; %position of centroid

% Plot fit with data.
centroid = centroid-(L/2);
plot(t/60,centroid);

I=get_last_monotonic(gradient(centroid),7);

ft = fittype('a*exp(-b*t)', 'independent', 't', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-L/2 0];
opts.StartPoint = [1 k2];
opts.Upper = [L/2 1];

% Fit model to data.
[fitresult, gof] = fit( t(I:end)'-t(I), centroid(I:end), ft, opts );
rate=coeffvalues(fitresult);
rate
hold on;

plot(t(I:end)/60, feval(fitresult,  t(I:end)'-t(I)));
hold off;
xlabel('Time (min)');
ylabel('Peak Position');
lgd = legend('Simulation','Fit', 'Location', 'north east');
%% Dispersion Profile

        %Jacobian
        fu=k2;
        fv=k3*C1^2;
        gu=-2*k2;
        gv=-k3*C1^2;
        
        nrange = [0:4*L];

        for i=1:length(nrange)
            %k=pi*n/L;
            k=pi*nrange(i)/L;
            k_2 = k.^2;
            AA=[fu-Da*k_2, fv;
            gu, gv-Db*k_2];
            lambda=eig(AA);
            [RE,I]=max(real(lambda));
            IM=imag(eig(I));
            A(i) = RE;
            A_im(i) = IM; 
        end
        
        figure(8);
        plot(nrange,A,'-o',nrange,A_im)
        title('Dispersion Relation','FontSize',12);
        xlabel('n')
        ylabel('growth rate')
        set(gca,'FontName', 'Helvetica','FontSize', 10);


%% Turing Space
% 
% figure(9)
% 
% cond=@(Da,Db,b,k2,k3) -Da*k3*b^2 +Db*k2 - 2*sqrt(Da*Db*(-k2*k3*b^2 + k3*b^2*2*k2));
% cond1 = @(b,k2,k3) k2/k3-b^2;
% cond2 = @(Da,Db,b,k2,k3) -k2*k3*b^2 + k3*b^2*2*k2;
% % Da = 0.01;
% % Db = 2;
% % b=1;
% 
% h = ezplot(@(k2,k3) cond(Da,Db,b,k2,k3));
% hold on
% h = ezplot(@(k2,k3) cond1(b,k2,k3));
% h = ezplot(@(k2,k3) cond2(Da,Db,b,k2,k3));
% set(h,'LineWidth',2);

% --------------------------------------------------------------
function [c,f,s] = pdes(~,~,u,DuDx)

c = [1; 1];
f = [Da;Db] .* DuDx; 
F = k3*u(2)*(u(1))^2;
s = [F-k2*u(1);-F+k4]; 
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



end