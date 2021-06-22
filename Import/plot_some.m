%plot data from multiple .mat files

clear all;

%%
%load files

z_1=load('Pal_0_ara.mat');
z_2=load('Pal_02_ara.mat');
z_3=load('Pal_dividing_30s.mat');
%z_4=load('nondiv2xbleach.mat');
%z_5=load('TolB-mCherry_nondividing_1s.mat');

%%
%violin plots

if size(z_1.D_median)==size(z_2.D_median)==size(z_3.D_median)

    figure(1)
    clf
    D_median=[z_1.D_median', z_2.D_median'; z_3.D_median'];%D_medians have to be of the same size
    violinplot(D_median);
    ylabel('Effective Diffusion constant (\mu m^2/s)')
    box on;
    set(gca, 'xtick', 1:6, 'xticklabels', {'z1','z2','z3','z4'});
    %ylim([0 0.02])

else
    disp('Cannot create violin plot, D_medians are of different size')
end
    

%%
%plot deff

figure(2)
clf
D=z_3.D;
CI=bootci(1000,@nanmedian,D');
shadederror(-1/2:0.02:1/2,nanmedian(D,2)',CI(1,:),CI(2,:),'WT')

hold on
D=z_1.D;
CI=bootci(1000,@nanmedian,D');
shadederror(-1/2:0.02:1/2,nanmedian(D,2)',CI(1,:),CI(2,:),'0% ara')


D=z_2.D;
CI=bootci(1000,@nanmedian,D');
shadederror(-1/2:0.02:1/2,nanmedian(D,2)',CI(1,:),CI(2,:),'0.2% ara')
ylim([0 3e-3])
%{
z_4.D_median(D_median>0.04)=NaN;%remove any bad cells
z_4.D=z_4.D(:,~isnan(z_4.D_median));
D=z_4.D;
CI=bootci(1000,@nanmedian,D');
shadederror(-1/2:0.02:1/2,nanmedian(D,2)',CI(1,:),CI(2,:),'Non-Dividing 2x')
D=z_5.D;
CI=bootci(1000,@nanmedian,D');
shadederror(0:0.02:1,nanmedian(D,2)',CI(1,:),CI(2,:),'Non-Dividing, 1s, 60s')
%}

hold off
legend;
xlabel('Relative Position')
ylabel('Effective diffusion constant (\mu m^2/s)')

%%
%plot kymogrpahs

guess=1e-3/(z_1.binfact*z_1.pixelsize)^2; %in units of binned pixels

data=z_1.avg;
[~,~,d2,~]=fitkymo(z_1.t,data,guess);

figure(101)
clf
subplot(2,1,1)
imagesc(-2:2:150,[-1/2,1/2],data)
xlabel('Time after bleaching (s)')
ylabel('Relative Position')
title('Averaged scaled data')

subplot(2,1,2)
sol=spatialFRAP(z_1.t,data(:,2),d2./data(:,1)/length(data(:,1)));
imagesc(-2:2:150,[-1/2,1/2],[data(:,1),sol])
xlabel('Time after bleaching (s)')
ylabel('Relative Position')
title('Fokker-Planck')


data=z_2.avg;
[~,~,d2,~]=fitkymo(z_2.t,data,guess);

figure(102)
clf
subplot(2,1,1)
imagesc(-2:2:150,[-1/2,1/2],data)
xlabel('Time after bleaching (s)')
ylabel('Relative Position')
title('Averaged scaled data')

subplot(2,1,2)
sol=spatialFRAP(z_2.t,data(:,2),d2./data(:,1)/length(data(:,1)));
imagesc(-2:2:150,[-1/2,1/2],[data(:,1),sol])
xlabel('Time after bleaching (s)')
ylabel('Relative Position')
title('Fokker-Planck')

%{
data=z_3.avg;
[~,~,d2,~]=fitkymo(z_3.t,data,guess);

figure(5)
clf
subplot(2,1,1)
imagesc(-5:5:100,[0,1],data)
xlabel('Time after bleaching (s)')
ylabel('Relative Position')
title('Averaged scaled data')

subplot(2,1,2)
sol=spatialFRAP(z_3.t,data(:,2),d2./data(:,1)/length(data(:,1)));
imagesc(-5:5:100,[0,1],[data(:,1),sol])
xlabel('Time after bleaching (s)')
ylabel('Relative Position')
title('Fokker-Planck')
%}
