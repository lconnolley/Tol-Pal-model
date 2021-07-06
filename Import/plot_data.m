function [diff] = plot_data(t,D_median,D,binfact,pixelsize,avg)

%violin plot of individual cells

x=-1/2:0.02:1/2;

figure(1)
clf;
D_median(D_median>1)=NaN;%remove any bad cells
D=D(:,~isnan(D_median));
violinplot(D_median');
ylabel('Effective Diffusion constant (\mu m^2/s)')
box on;

%%
%average effective diffusion coefficient of individual cells

figure(2)
clf
CI=bootci(1000,@nanmedian,D');
shadederror(x,nanmedian(D,2)',CI(1,:),CI(2,:),'data')
xlabel('Relative Position')
ylabel('Effective diffusion constant (\mu m^2/s)')
%ylim([0 14*10^(-4)])


%%
%kymographs of average cell

guess=1e-2/(binfact*pixelsize)^2; %in units of binned pixels

data=avg;
[~,~,d2,~]=fitkymo(t,data,guess);
%d1*(binfact*pixelsize)^2; %Fickian
diff=d2*(pixelsize)^2; %Fokker Planck

%set colourbar to be the same across kymogrpahs
sol=spatialFRAP(t,data(:,2),d2./data(:,1)/length(data(:,1)));
bottom = min(min(min(data)),min(min(sol)));
top  = max(max(max(data)),max(max(sol)));

delta=t(2)-t(1);
t1=[-delta; t']';

figure(3)
clf
subplot(2,1,1)
imagesc(t1,[-1/2,1/2],data)
title('Average scaled data')
caxis manual
caxis([bottom top])

subplot(2,1,2)
imagesc(t1,[-1/2,1/2],[data(:,1),sol])
title('Fokker-Planck')
caxis manual
caxis([bottom top]);

%%
%average diffusion coefficient

%{
figure(5)
clf
plot(x,diff./data(:,1)/length(data(:,1)))
%}

end