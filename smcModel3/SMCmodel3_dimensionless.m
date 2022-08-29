function SMCmodel3_dimensionless

m = 0;
L = 1;
x = 0:0.01:L;
t = 0:0.01:2;
Du = 0.3;
Dv = 0.012;


a = 3.75;
b = 0.0039;
Gamma = 1200;

d = Du/Dv;


u0 = (b+1)/(a+b+1);
v0 = 1-u0;


rr1=randn(size(x));
rr2=randn(size(x));
rr2=sum(rr1)*rr1/sum(rr2);
IC = [u0*(1+0.1*rr1);v0*(1-0.1*rr2)];
% IC = [rr1;rr2];


options=odeset('RelTol',1e-8,'AbsTol',1e-16);
sol = pbcpdeSolver(@fpde,@ic,x,t);
u1 = sol(:,:,1);
u2 = sol(:,:,2);



%%Plot
figure(1);
clf;
subplot(4,1,1);
plot(x,u2(end,:),x,u1(end,:))
title('Solution at t = end');
xlabel('Distance x');
%ylim([0,2*C]);

for n=0:40
an(n+1)=2*trapz(x,u2(end,:).*cos(pi*n*x/L)/1000)/L;
end
an = an(2:16);
[val,idx] = max(abs(an))

figure(2); 
clf;
subplot(2,1,1);
imagesc(t/60,x,u1');
%set(p1,'LineStyle','none');
title('u(x,t)')
ylabel('Position (micro m)')
xlabel('Time (min)')
%zlim([0,5]);

subplot(2,1,2);
imagesc(t/60,x,u2');
%set(p2,'LineStyle','none')
title('v(x,t)')
ylabel('Position (micro m)')
xlabel('Time (min)')
%zlim([0,15]);



% Centroid rate
figure(3)

centroid=diag(1./sum(u2,2))*u2*x'; %position of centroid

% Plot fit with data.
centroid = centroid-(L/2);
plot(t/60,centroid);

% I=get_last_monotonic(gradient(centroid),7);
I=10;
% Centroid rate fit
ft = fittype('a*exp(-b*t)', 'independent', 't', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-L/2 0];
opts.StartPoint = [1 0.2];
opts.Upper = [L/2 1];

% Fit model to data.
[fitresult, gof] = fit( t(I:end)'-t(I), centroid(I:end), ft, opts );
rate=coeffvalues(fitresult);
rate
hold on;
plot(t(I:end)/60, feval(fitresult,  t(I:end)'-t(I)))
hold off;
xlabel('time (min)');
ylabel('centroid')
% lgd = legend('Centroid Position','Exponential Fit', 'Location', 'north east');

%% Dispersion Profile

        %Jacobian
        fu = -a-2*a*u0-b;
        fv = 1-2*a*u0;
        gu = a+2*a*u0;
        gv = 2*a*u0-1-b;
        
        nrange = [0:8*L];

        for i=1:length(nrange)
            %k=pi*n/L;
            k=pi*nrange(i);
            k_2 = k.^2;
            AA =  Gamma*[fu, fv;
            gu, gv] - [d*k_2, 0;0, k_2];
            lambda=eig(AA);
            [RE,I]=max(real(lambda));
            IM=imag(eig(I));
            A(i) = RE;
            A_im(i) = IM; 
        end
        pred = find(A==max(A))-1
        figure(9);
        plot(nrange,A,'-o',nrange,A_im)
        title('Dispersion Relation','FontSize',12);
        xlabel('n')
        ylabel('growth rate')
        set(gca,'FontName', 'Helvetica','FontSize', 10);
%--------------------------------------------------------------
function [c,f,s] = fpde(~,~,u,DuDx)

c = [1; 1];
f = [d;1] .*DuDx; 
F = a*u(1)*(u(1)+u(2))^2-u(2);%
s = [Gamma*(-F + b*(1-u(1))); Gamma*(F-b*u(2))]; 
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