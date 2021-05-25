clear all;
name='div2xbleach.mat';
load(name);
name='Non-dividing';

%%
x=-1/2:0.02:1/2;

%individual cells
figure(1)
clf;
%D_median(D_median>0.1)=NaN;%remove any bad cells
D=D(:,~isnan(D_median));
violinplot([D_median']);
ylabel('Effective Diffusion constant (\mu m^2/s)')
box on;

figure(2)
clf
CI=bootci(1000,@nanmedian,D');
shadederror(x,nanmedian(D,2)',CI(1,:),CI(2,:),name)
legend;
xlabel('Relative Position')
ylabel('Effective diffusion constant (\mu m^2/s)')
%ylim([0 14*10^(-4)])

%figure(3)
%clf
%plot(0:0.02:1,1./data(:,1)/length(data(:,1)))

%%
%averaged cell

guess=1e-2/(binfact*pixelsize)^2; %in units of binned pixels

data=avg;
[d1,fval1,d2,fval2]=fitkymo(t,data,guess);
%d1*(binfact*pixelsize)^2 %Fickian
%fval1
diff=d2*(adj_pixelsize*pixelsize)^2 %Fokker Planck
fval2

%colourbar
sol=spatialFRAP(t,data(:,2),d2./data(:,1)/length(data(:,1)));
bottom = min(min(min(data)),min(min(sol)));
top  = max(max(max(data)),max(max(sol)));

t1=[-2; t']';

figure(4)
clf
subplot(2,1,1)
imagesc(t1,[-1/2,1/2],data)
title('Averaged scaled data')

subplot(2,1,2)
imagesc(t1,[-1/2,1/2],[data(:,1),sol])
title('Fokker-Planck')

caxis manual
caxis([bottom top]);

%average diffusion coefficient
figure(5)
clf
plot(x,diff./data(:,1)/length(data(:,1)))
