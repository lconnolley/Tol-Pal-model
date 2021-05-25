clear all;

z_1=load('div3xbleach.mat');
z_2=load('div2xbleach.mat');
z_3=load('nondiv3xbleach.mat');
z_4=load('nondiv2xbleach.mat');
%z_5=load('TolB-mCherry_nondividing_1s.mat');

%%

figure(1)
clf;
D_median=z_1.D_median';
violinplot(D_median); 
ylabel('Effective Diffusion constant (\mu m^2/s)')
box on;
set(gca, 'xtick', 1:6, 'xticklabels', {'TolB Dividing'});
%ylim([0 0.02])

figure(2)
clf;
D_median=z_2.D_median';
violinplot(D_median);
ylabel('Effective Diffusion constant (\mu m^2/s)')
box on;
set(gca, 'xtick', 1:6, 'xticklabels', {'TolB non-dividing'});
%ylim([0 0.02])

%{
figure(3)
clf;
D_median=z_2.D_median';
violinplot(D_median);
ylabel('Effective Diffusion constant (\mu m^2/s)')
box on;
set(gca, 'xtick', 1:6, 'xticklabels', {'Non-Dividing 2s, 150s'});
ylim([0 0.02])

figure(4)
clf
D_median=z_4.D_median';
violinplot(D_median)
ylabel('Effective Diffusion constant (\mu m^2/s)')
box on;
set(gca, 'xtick', 1:6, 'xticklabels', {'Non-Dividing 2s 20s, 300s'});
ylim([0 0.02])
%}

figure(5)
clf
D=z_1.D;
CI=bootci(1000,@nanmedian,D');
shadederror(-1/2:0.02:1/2,nanmedian(D,2)',CI(1,:),CI(2,:),'Dividing 3x')

hold on
D=z_2.D;
CI=bootci(1000,@nanmedian,D');
shadederror(-1/2:0.02:1/2,nanmedian(D,2)',CI(1,:),CI(2,:),'Dividing 2x')

D=z_3.D;
CI=bootci(1000,@nanmedian,D');
shadederror(-1/2:0.02:1/2,nanmedian(D,2)',CI(1,:),CI(2,:),'Non-Dividing 3x')

z_4.D_median(D_median>0.04)=NaN;%remove any bad cells
z_4.D=z_4.D(:,~isnan(z_4.D_median));
D=z_4.D;
CI=bootci(1000,@nanmedian,D');
shadederror(-1/2:0.02:1/2,nanmedian(D,2)',CI(1,:),CI(2,:),'Non-Dividing 2x')

ylim([0 8e-4])

%{
D=z_5.D;
CI=bootci(1000,@nanmedian,D');
shadederror(0:0.02:1,nanmedian(D,2)',CI(1,:),CI(2,:),'Non-Dividing, 1s, 60s')
%}
legend;
xlabel('Relative Position')
ylabel('Effective diffusion constant (\mu m^2/s)')

%%

guess=1e-3/(z_1.binfact*z_1.pixelsize)^2; %in units of binned pixels

data=z_1.avg;
[d1,fval1,d2,fval2]=fitkymo(z_1.t,data,guess);
%F=d1*(z_1.binfact*z_1.pixelsize)^2 %Fickian
%fval1
FP=d2*(z_1.binfact*z_1.pixelsize)^2 %Fokker Planck
fval2

figure(6)
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
[d1,fval1,d2,fval2]=fitkymo(z_2.t,data,guess);
%F=d1*(z_2.binfact*z_2.pixelsize)^2 %Fickian
%fval1
FP=d2*(z_2.binfact*z_2.pixelsize)^2 %Fokker Planck
fval2

figure(7)
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
[d1,fval1,d2,fval2]=fitkymo(z_3.t,data,guess);
%F=d1*(z_3.binfact*z_3.pixelsize)^2 %Fickian
%fval1
FP=d2*(z_3.binfact*z_3.pixelsize)^2 %Fokker Planck
fval2

figure(8)
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
