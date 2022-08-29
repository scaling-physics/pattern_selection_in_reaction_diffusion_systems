%%SMCmodel3_two_peaks_out
% addpath('/Users/srik/Documents/MATLAB/PlotPub-master/lib');

s=[];
load('two_peak_numeric_traj_1.mat');
for shift = 160
   for shift1 = 160
[peak1,peak2,utotal] = SMCmodel3_two_peaks_utotal(shift,shift1);

        if(length(peak1) == length(utotal))
            s1.mat = [peak1, peak2, utotal];
            s = [s,s1];
        end
save('two_peak_numeric_traj_2.mat','s');
    end
end

%% Plot

% load('two_peak_numeric_traj_2.mat');
load('two_inverse_sinks_mass.mat');
sa=[];

for i=1:length(s)
    
    if abs(s(i).mat(end,1)+s(i).mat(end,2)) < 0.01 && s(i).mat(end,3) < 0.188 
        sa = [sa; s(i).mat];
    end
end

set(0,'defaultAxesFontSize',12)
figure(1)
x = sa(:,4);
y = sa(:,5);
z = sa(:,3);

% z = normalize(z,'range');
% x = [x;-x];
% y = [y;-y];
% z = [z;z];




% F = scatteredInterpolant(x,y,z);
% F.Method = 'linear';


[xq,yq] = meshgrid(0:-0.001:-0.5, 0:0.001:0.5);
% [xq,yq] = meshgrid(min(x):0.001:max(x), min(y):0.001:max(y));
% [xq,yq] = meshgrid(-0.42:0.001:-0.1, 0.1:0.001:0.42);
% cq = F(xq,yq);
% h = pcolor(xq,yq,cq);
% h.EdgeColor = 'none';
% c = ;
% colormap(parula)
% zq = interp3(x,y,z,xq,yq);
zq = griddata(x,y,z,xq,yq,'linear');
[c,h] = contourf(xq,yq,zq,'LineWidth',2);
% clabel(c,h);
% c.Label.Interpreter = 'latex';
% c.Label.FontSize = 24;
cb = colorbar();
cb.Label.String = 'M/c';
xlabel('x_1/L');
ylabel('x_2/L');
hold on
for i=1:2:length(s)
    
%     if abs(abs(s(i).mat(end,1)) - 0.25) < 0.01 && s(i).mat(1,3) < 0.19
        plot(s(i).mat(:,1), s(i).mat(:,2),'w','LineWidth',2);

%     end
end
% load('two_sinks_mass.mat')
for i=1:2:length(s)
    
%     if abs(abs(s(i).mat(end,1)) - 0.25) < 0.01 && s(i).mat(1,3) < 0.19
        plot(s(i).mat(:,3), s(i).mat(:,4),'r','Linestyle','--','LineWidth',2);

%     end
end
hold off
xlim([-0.5 0])
ylim([0 0.5])

set(gca, 'xdir', 'reverse','Xtick',-.5:0.25:0,'Ytick',0:0.25:0.5)

box on


% export_fig('/Users/srik/Documents/Pattern_selection_images/Figure 4/utotal_two_peaks_numeric_2D_Dv*0.1.pdf');

% export_fig('/Users/srik/Documents/Pattern_selection_images/Figure 4/utotal_two_peaks_anayltic_2D_Dv*0.1.pdf');