%fit to Pal dividing and non-dividing data
clear all

%%
%load experimental data

P_d=load('../Import/Pal_dividing_30s.mat');
P_nd=load('../Import/Pal_nondividing_30s.mat');

%% preprocess experimental data so in correct form to fit to

%renormalise data to tolb/pal concentration after bleach, change so 
%normalised to 1 on x, multiply by average concentration of total tolb/pal 
%in a cell (uM), multiply by area under bleached curve
div=(P_d.avg/trapz(P_d.avg(:,1)))*99.6e3*P_d.factor;
non_div=(P_nd.avg/trapz(P_nd.avg(:,1)))*99.6e3*P_nd.factor;

%calculate average cell length in um
lngth=cellfun('length',P_d.cells);
Pd_lngth=median(lngth)/50;
lngth=cellfun('length',P_nd.cells);
Pnd_lngth=median(lngth)/50;

%calculate average effective diffusion coefficients
Pnd_data=P_nd.avg;
[d2]=fitkymo(P_nd.t,Pnd_data,0.01);
Pnd_diff=d2*(Pnd_lngth*P_nd.pixelsize)^2;
Pd_data=P_d.avg;
[d2]=fitkymo(P_d.t,Pd_data,0.01);
Pd_diff=d2*(Pd_lngth*P_d.pixelsize)^2;

deff_d=Pd_diff./Pd_data(:,1)/length(Pd_data(:,1));
deff_d1=[0.000788*ones(9,1); deff_d(10:42,:); 0.0008582*ones(9,1)]';%remove weird edges from dividing deff
deff_nd=Pnd_diff./Pnd_data(:,1)/length(Pnd_data(:,1));

%% fit model to experimental data

%guess=[a, b, beta0]
guess=[0.05 2.7 20]; %a=Dc-Db, b=Dc/Db, beta0

%tic
[d,fval]=fitkymo_pal(div,non_div,deff_d1,deff_nd,guess);
Dc=(d(2)*d(1))/(d(2)-1);
Db=d(1)/(d(2)-1);
beta0=d(3);
%toc

save('fit_parameter_new.mat','d','Dc','Db');%save parameter values found

%% find prebleach steady state solution

[w1,w2,w3,w4]=steady_state_d(d(1),d(2),d(3),Pd_lngth);
[w10,w20,w30,w40]=steady_state_nd(d(1),d(2),d(3),Pnd_lngth);

%% plot pal dividing and non-dividing kymographs
%N.B. these are not normalised properly

figure(1)
clf
subplot(2,1,1) 
imagesc(-30:30:600,[-1/2,1/2],div)
title('Average scaled data - Pal dividing')
subplot(2,1,2)
f=spatialFRAP_Pal_d(d(1),d(2),d(3));
imagesc(-30:30:600,[-1/2,1/2],f)
title('Computational - Pal dividing')

trapz(f);

figure(2)
clf
subplot(2,1,1)
imagesc(-30:30:600,[-1/2,1/2],non_div)
title('Average scaled data - Pal non-dividing')
subplot(2,1,2)
g=spatialFRAP_Pal_nd(d(1),d(2),d(3));
imagesc(-30:30:600,[-1/2,1/2],g)
title('Computational - Pal non-dividing')

trapz(g);



