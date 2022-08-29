clear;
filename = 'SMCmodel3_turing_space_periodic_numerical.mat';
z=load(filename);

x=unique(z.mat(:,1));
y=unique(z.mat(:,2));

[X,Y]=ndgrid(x,y);

runs=cell(size(X));
for i=1:length(z.mat(:,1))
    x0=z.mat(i,1);
    y0=z.mat(i,2);
    Ix=find(abs(x-x0)<0.001);
    Iy=find(abs(y-y0)<0.01);
    runs{Ix,Iy}=z.mat(i,3:end);
end

modes=NaN(size(runs));
for i=1:size(runs,1)
    for j=1:size(runs,2)
        if ~isempty(runs{i,j})
            modes(i,j)=mode(runs{i,j});
        end
    end
end


figure(1)

h=imagesc(x,y,modes');
set(gca, 'YDir','normal')
colorbar;
set(h,'AlphaData',~isnan(modes'));

