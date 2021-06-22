clear all

%%
%load in SpatialFRAP files

Pd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Fitting/Pal_dividing.mat');
Pnd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/non-dividing.mat');
A=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolA_dividing.mat');
B=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolB_dividing.mat');

%Bd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_dividing_5s.mat');%peaks and no peaks (no selection)
Bd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_dividing_nopeaks.mat');%no peaks
Bnd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_nondiv_2s_150s.mat');

%%
%calculate average length of cells in um

no_segments=50;

lngth=cellfun('length',Pd.cells);
Pd_lngth=median(lngth);
Pd_adjust=Pd_lngth/no_segments;

lngth=cellfun('length',Pnd.cells);
Pnd_lngth=median(lngth);
Pnd_adjust=Pnd_lngth/no_segments;

lngth=cellfun('length',A.cells);
A_lngth=median(lngth);
A_adjust=A_lngth/no_segments;

lngth=cellfun('length',B.cells);
B_lngth=median(lngth);
B_adjust=B_lngth/no_segments;

lngth=cellfun('length',Bd.cells);
Bd_lngth=median(lngth);
Bd_adjust=Bd_lngth/no_segments;

lngth=cellfun('length',Bnd.cells);
Bnd_lngth=median(lngth);
Bnd_adjust=Bnd_lngth/no_segments;

%%
%fit to experimental SpatialFRAP to find deff for Pal data

x=-1/2:0.02:1/2;
guess=0.001;

Pnd_data=Pnd.avg;
[~,~,d2,~]=fitkymo_avg(Pnd.t,Pnd_data,guess);
Pnd_diff=d2*(Pnd_adjust*Pnd.pixelsize)^2;

Pd_data=Pd.avg;
[~,~,d2,~]=fitkymo_avg(Pd.t,Pd_data,guess);
Pd_diff=d2*(Pd_adjust*Pd.pixelsize)^2;

A_data=A.avg;
[~,~,d2,~]=fitkymo_avg(A.t,A_data,guess);
A_diff=d2*(A_adjust*A.pixelsize)^2;

B_data=B.avg;
[~,~,d2,~]=fitkymo_avg(B.t,B_data,guess);
B_diff=d2*(B_adjust*B.pixelsize)^2;

%%
%plot effective diffusion coefficient for average of cells for Pal data

figure(1)
clf 
plot(x,Pnd_diff./Pnd_data(:,1)/length(Pnd_data(:,1)),'DisplayName','Non-dividing','LineWidth',2)
hold on
plot(x,Pd_diff./Pd_data(:,1)/length(Pd_data(:,1)),'DisplayName','Dividing','LineWidth',2)
plot(x,A_diff./A_data(:,1)/length(A_data(:,1)),'DisplayName','tolA','LineWidth',2)
plot(x,B_diff./B_data(:,1)/length(B_data(:,1)),'DisplayName','tolB','LineWidth',2)
%hold off
legend
%ylim([0 1e-3])
title('Average')
xlabel('Relative position')
ylabel('Effective diffusion constant (\mu m^2/s)')
%%
%fit to experimental SpatialFRAP to find deff for TolB data

Bd_data=Bd.avg;
[~,~,d2,~]=fitkymo_avg(Bd.t,Bd_data,guess);
Bd_diff=d2*(Bd_adjust*Bd.pixelsize)^2;

Bnd_data=Bnd.avg;
[~,~,d2,~]=fitkymo_avg(Bnd.t,Bnd_data,guess);
Bnd_diff=d2*(Bnd_adjust*Bnd.pixelsize)^2;

%%
%plot effective diffusion coefficient for average of cells for TolB data

figure(2)
clf
plot(x,Bnd_diff./Bnd_data(:,1)/length(Bnd_data(:,1)),'DisplayName','Non-dividing','LineWidth',2)
hold on
plot(x,Bd_diff./Bd_data(:,1)/length(Bd_data(:,1)),'DisplayName','Dividing','LineWidth',2)
hold off
legend
%ylim([0 0.07])
%title('TolB average')
xlabel('Relative position')
ylabel('Effective diffusion constant (\mu m^2/s)')

