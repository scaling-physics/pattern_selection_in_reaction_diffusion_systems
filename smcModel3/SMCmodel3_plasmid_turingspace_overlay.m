L = 4;
m = 0;

gamma=3.6;
s = 1;
Dv = 0.012*s;
% Dv = 0.0001;
D = 0.3;

c = 300;

%% Plasmid parameters mapping

n = 1:100;
%% Sweep
d = D/Dv;
mat = zeros(1,3);    
cond=@(d,b,a) d*(b+1).*(a-b-1)./(a+b+1)+(-(a.^2+b.^2+4*a.*b+3.*a+b)./(a+b+1))-2*sqrt(d*(b).*(a+b+1));
% cond = @(b,a) 4*(b+1)-a;
% cond1 = @(d,b,a) (1-a*0.5*(1+sqrt(1-4*(b+1)/a)) + d*(b+1)).^2 - 4*a*b*(0.5*(1+sqrt(1-4*(b+1)/a)) - 2*(b+1)/a);

% h = ezplot(@(b,a) cond(d,b,a),[0,10,0,5000]);

for j = 1:20:200000
    for m = 0:2:250
       
        
       k = 0.001*j;


       b = k^2/L^2*D/gamma;
       delta = b*gamma;
       a = 0.1*m;
%        a = 5;
       if cond(d,b,a) > 0
            
            
            
%             delta = gamma*b;
            beta = a*gamma/c^2;
%             k = L*sqrt(delta/D);
            gp = coth(k./(2.*n));
            rho = 6.*sqrt(Dv).*delta.*sqrt(delta + gamma)./beta./L;
            muprime1 = c.*delta./(gp.*k);
            muprimep = muprime1 + sqrt(muprime1.^2 - (2.*delta)./(gp.*k).*rho);
            muprimem = muprime1 - sqrt(muprime1.^2 - (2.*delta)./(gp.*k).*rho);
            
            M1 = c - n/(delta).*muprimep;

            M2 = c - n/(delta).*muprimem;
            
            A = 1./M1(~imag(M1));
           
            nmax = find(A == max(A));

            mat = [mat; k a nmax];
       
        end
    end
end
%%
% filename = "/Users/srik/Documents/MATLAB/flow/data/raw_data/SMCmodel2_turing_spaceoverlay_plasmid_Dv0.5.mat";
% save(filename);

%% Plot

figure(1)

scatter(mat(:,1),mat(:,2),10,mat(:,3),'filled');
p1=xlabel('$\kappa$');
p=ylabel('$a$');
set(p,'Interpreter','latex','FontSize',16)
set(p1,'Interpreter','latex','FontSize',16);
cb = colorbar();
% ylabel(cb, 'Dominant unstable mode')
ylabel(cb, 'Number of peaks in dominant mode')
% colorbar('XTickLabel',{'0','1','2','3','4','5'})
numcolors = max(mat(:,3))+1;
colormap(parula(numcolors));
caxis([-0.5 numcolors-0.5]);


% figure(2)
% scatter(mat(:,1),mat(:,3),30,'filled')
% ylim([0 max(mat(:,3))+2])