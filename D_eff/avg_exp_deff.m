clear all

Pd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Fitting/Pal_dividing.mat');
Pnd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/non-dividing.mat');
A=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolA_dividing.mat');
B=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolB_dividing.mat');

Bd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_dividing_nopeaks.mat');
Bnd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_nondiv_2s_150s.mat');

for i=1:length(Pd.cells)
    data=Pd.cells{i};
    data=data(3:end-2,:);
    npixels(i)=length(data(:,1));
end
Pd_lngth=median(npixels);
Pd_adjust=Pd_lngth/51;

for i=1:length(Pnd.cells)
    data=Pnd.cells{i};
    data=data(3:end-2,:);
    npixels(i)=length(data(:,1));
end
Pnd_lngth=median(npixels);
Pnd_adjust=Pnd_lngth/51;

for i=1:length(A.cells)
    data=A.cells{i};
    data=data(3:end-2,:);
    npixels(i)=length(data(:,1));
end
A_lngth=median(npixels);
A_adjust=A_lngth/51;

for i=1:length(B.cells)
    data=B.cells{i};
    data=data(3:end-2,:);
    npixels(i)=length(data(:,1));
end
B_lngth=median(npixels);
B_adjust=B_lngth/51;

for i=1:length(Bd.cells)
    data=Bd.cells{i};
    data=data(3:end-2,:);
    npixels(i)=length(data(:,1));
end
Bd_lngth=median(npixels);
Bd_adjust=Bd_lngth/51;

for i=1:length(Bnd.cells)
    data=Bnd.cells{i};
    data=data(3:end-2,:);
    npixels(i)=length(data(:,1));
end
Bnd_lngth=median(npixels);
Bnd_adjust=Bnd_lngth/51;


x=-1/2:0.02:1/2;
guess=1e-2/(Pd_adjust*Pd.pixelsize)^2;

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


figure(1)
clf
plot(x,Pnd_diff./Pnd_data(:,1)/length(Pnd_data(:,1)),'DisplayName','Non-dividing','LineWidth',2)
hold on
plot(x,Pd_diff./Pd_data(:,1)/length(Pd_data(:,1)),'DisplayName','Dividing','LineWidth',2)
plot(x,A_diff./A_data(:,1)/length(A_data(:,1)),'DisplayName','tolA','LineWidth',2)
plot(x,B_diff./B_data(:,1)/length(B_data(:,1)),'DisplayName','tolB','LineWidth',2)
legend
title('Average')
xlabel('Relative Position')
ylabel('Effective diffusion constant (\mu m^2/s)')

figure(2)
clf
plot(x,median(Pnd.D,2),'DisplayName','Non-dividing','LineWidth',2)
hold on
plot(x,median(Pd.D,2),'DisplayName','Dividing','LineWidth',2)
plot(x,median(A.D,2),'DisplayName','tolA','LineWidth',2)
plot(x,median(B.D,2),'DisplayName','tolB','LineWidth',2)
legend
title('Individual')
