function brusselator_dimensionless

m = 0;
L= 1;
x = 0:0.01:L;
t = 0:1:100;
Du = 1;
Dv = 0.01;



a = 0.5;
Gamma = 1500;
d = Du/Dv;

b = 0;
u0 = 1/a;
v0 = 1;


rr1=randn(size(x));
rr2=randn(size(x));
rr2=sum(rr1)*rr1/sum(rr2);
IC = [u0*(1+0.01*rr1);v0*(1-0.01*rr2)];


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
% 
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
        fu = -a;
        fv = -1;
        gu = a;
        gv = 1-b;
        
        nrange = [0:30*L];

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
        figure(8);
        plot(nrange,A,'-o',nrange,A_im)
        title('Dispersion Relation','FontSize',12);
        xlabel('n')
        ylabel('growth rate')
        set(gca,'FontName', 'Helvetica','FontSize', 10);
%--------------------------------------------------------------
function [c,f,s] = pdes(~,~,u,DuDx)

c = [1; 1];
f = [d;1] .* DuDx; 
F = u(1)*(u(2))^2;
s = [Gamma*(-a*F+u(2));Gamma*(a*F-u(2)+b*(1-u(2)))]; 
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