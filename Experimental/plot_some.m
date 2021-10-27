%plot data from multiple .mat files, can create plots for figure 1, figure 
%4 - S1, figure 5(e), figure 6

clear all

%% load files

%z_1 = load('../Import/TolB_nondiv_2s.mat');
%z_2 = load('../Import/TolB_dividing_nopeaks.mat');
z_1 = load('../Import/Pal_nondividing_30s.mat');
z_2 = load('../Import/Pal_dividing_30s.mat');
z_3 = load('../Import/tolA_dividing.mat');
z_4 = load('../Import/tolB_dividing.mat');
%z_1=load('../Import/Pal_0_ara.mat');
%z_2=load('../Import/Pal_02_ara.mat');
%z_5=load('../Import/TolB-mCherry_nondividing_1s.mat');

%% violin plots

if size(z_1.D_median)==size(z_2.D_median)

    figure(1)
    clf
    %z_1.D_median(z_1.D_median>0.0006)=NaN;
    %z_2.D_median(z_2.D_median>0.0006)=NaN;
    D_median=[z_1.D_median', z_2.D_median'];%; z_4.D_median'];%D_medians have to be of the same size
    violinplot(D_median);
    ylabel('Effective Diffusion constant (\mu m^2/s)')
    box on;
    set(gca, 'xtick', 1:6, 'xticklabels', {'z1','z2','z3','z4'});
    %ylim([0 0.02])

else
    disp('Cannot create violin plot, D_medians are of different size')
end
    

%% plot deff

figure(2)
clf
D=z_1.D;
CI=bootci(1000,@nanmedian,D');
shadederror(-1/2:0.02:1/2,nanmedian(D,2)',CI(1,:),CI(2,:),'Non-dividing')

hold on
D=z_2.D;
CI=bootci(1000,@nanmedian,D');
shadederror(-1/2:0.02:1/2,nanmedian(D,2)',CI(1,:),CI(2,:),'Dividing')


D=z_3.D;
CI=bootci(1000,@nanmedian,D');
shadederror(-1/2:0.02:1/2,nanmedian(D,2)',CI(1,:),CI(2,:),'tolA')

%D=z_4.D;
%CI=bootci(1000,@nanmedian,D');
%shadederror(-1/2:0.02:1/2,nanmedian(D,2)',CI(1,:),CI(2,:),'tolB')
%ylim([0 2e-3])

ylim([0 2e-3])
hold off
legend;
xlabel('Relative Position')
ylabel('Effective diffusion constant (\mu m^2/s)')

%% plot kymogrpahs

guess=1e-3/(z_1.binfact*z_1.pixelsize)^2; %in units of binned pixels

data=z_1.avg;
[~,~,d1,~]=fitkymo(z_1.t,data,guess);

data2=z_2.avg;
[~,~,d2,~]=fitkymo(z_2.t,data2,guess);

figure(101)
clf
subplot(2,1,1)
imagesc(-2:2:150,[-1/2,1/2],data)
xlabel('Time after bleaching (s)')
ylabel('Relative Position')
title('Averaged scaled data')
colorbar
subplot(2,1,2)
sol=spatialFRAP(z_1.t,data(:,2),d1./data(:,1)/length(data(:,1)));
imagesc(-2:2:150,[-1/2,1/2],[data(:,1),sol])
xlabel('Time after bleaching (s)')
ylabel('Relative Position')
title('Fokker-Planck')
colorbar

figure(102)
clf
subplot(2,1,1)
imagesc(-2:2:150,[-1/2,1/2],data2)
xlabel('Time after bleaching (s)')
ylabel('Relative Position')
title('Averaged scaled data')
colorbar
subplot(2,1,2)
sol=spatialFRAP(z_2.t,data2(:,2),d2./data2(:,1)/length(data2(:,1)));
imagesc(-2:2:150,[-1/2,1/2],[data2(:,1),sol])
xlabel('Time after bleaching (s)')
ylabel('Relative Position')
title('Fokker-Planck')
data2=z_2.avg;
[~,~,d2,~]=fitkymo(z_2.t,data2,guess);
colorbar

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
