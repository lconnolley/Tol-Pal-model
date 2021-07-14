%plot experimental and model kymographs for pal, TolB, and mutants
clear

%%
%load experimental data

B_d=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_dividing_nopeaks.mat');
P_d=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/Pal_dividing_30s.mat');%Pal_dividing.mat');

B_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_nondiv_2s.mat');
P_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/Pal_nondividing_30s.mat');%/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/non-dividing.mat');

A=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolA_dividing.mat');
B=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolB_dividing.mat');

%%
%parameters found from fitting
z=load('fit_parameter.mat');
%d = [a, b, beta0], a=Dc-Db, b=Dc/Db, beta0 
d=z.d;

xs=-1/2:0.005:1/2;
%%
%preprocess experimental data into same units as model data
xe=-1/2:0.02:1/2;

%renormalise data to tolb/pal concentration after bleach, change so 
%normalised to 1 on x, multiply by average concentration of total tolb/pal 
%in a cell (uM), multiply by area under bleached curve
tolb_d=(B_d.avg/trapz(xe,B_d.avg(:,1)));%*9.96e5*B_d.factor;                 %not sure why e5 not e3
pal_d=(P_d.avg/trapz(xe,P_d.avg(:,1)));%*99.6e5*P_d.factor;
tolb_nd=(B_nd.avg/trapz(xe,B_nd.avg(:,1)));%*9.96e5*B_nd.factor;
pal_nd=(P_nd.avg/trapz(xe,P_nd.avg(:,1)));%*99.6e5*P_nd.factor;
tolA=(A.avg/trapz(xe,A.avg(:,1)));%*99.6e5*A.factor;   %beta0=0
tolB=(B.avg/trapz(xe,B.avg(:,1)));%*99.6e5*B.factor;   %alpha=0

%%
%TolB dividing kymographs for experimental and model

f1=spatialFRAP_TolB_d(d(1),d(2),d(3));
f1=f1./trapz(xs,f1(:,1));

%bottom = min(min(min(tolb_d.*factor)),min(min(f1)));
%top  = max(max(max(tolb_d.*factor)),max(max(f1)));

figure(1)
clf
subplot(2,1,1)
imagesc(-2:2:150,[-1/2,1/2],tolb_d)
title('Average scaled data - TolB (dividing)')
%caxis manual
%caxis([bottom top]);
colorbar
subplot(2,1,2)
imagesc(-2:2:150,[-1/2,1/2],f1)
title('Computational - TolB (dividing)')
%caxis manual
%caxis([bottom top]);
colorbar

trapz(xe,tolb_d);
trapz(xs,f1);

%%
%Pal dividing kymographs 

g1=spatialFRAP_Pal_d(d(1),d(2),d(3));
g1=g1./trapz(xs,g1(:,1));

bottom = min(min(min(pal_d)),min(min(g1)));
top  = max(max(max(pal_d)),max(max(g1)));

figure(2)
clf
subplot(2,1,1)
imagesc(-30:30:600,[-1/2,1/2],pal_d)
title('Average scaled data - Pal (dividing)')
caxis manual
caxis([bottom top]);
colorbar
subplot(2,1,2)
imagesc(-30:30:600,[-1/2,1/2],g1)
title('Computational - Pal (dividing)')
caxis manual
caxis([bottom top]);
colorbar

trapz(xs,g1);
trapz(xe,pal_d);

%%
%TolB non-dividing kymographs

f2=spatialFRAP_TolB_nd(d(1),d(2),d(3));
f2=f2./trapz(xs,f2(:,1));

bottom = min(min(min(tolb_nd)),min(min(f2)));
top  = max(max(max(tolb_nd)),max(max(f2)));

figure(3)
clf
subplot(2,1,1)
imagesc(-2:2:150,[-1/2,1/2],tolb_nd)
title('Average scaled data - TolB (non-dividing)')
caxis manual 
caxis ([bottom top])
colorbar
subplot(2,1,2)
imagesc(-2:2:150,[-1/2,1/2],f2)
title('Computational - TolB (non-dividing)')
caxis manual
caxis ([bottom top])
colorbar

%%
%Pal non-dividing kymographs

g2=spatialFRAP_Pal_nd(d(1),d(2),d(3));
g2=g2./trapz(xs,g2(:,1));

bottom = min(min(min(pal_nd)),min(min(g2)));
top  = max(max(max(pal_nd)),max(max(g2)));

figure(4)
clf
subplot(2,1,1)
imagesc(-30:30:600,[-1/2,1/2],pal_nd)
title('Average scaled data - Pal (non-dividing)')
caxis manual
caxis ([bottom top])
colorbar
subplot(2,1,2)
imagesc(-30:30:600,[-1/2,1/2],g2)
title('Computational - Pal (non-dividing)')
caxis manual
caxis ([bottom top])
colorbar

%%
%tolA kymographs

f3=spatialFRAP_tolA_d(d(1),d(2),0);
f3=f3./trapz(xs,f3(:,1));

bottom = min(min(min(tolA)),min(min(f3)));
top  = max(max(max(tolA)),max(max(f3)));

figure(5)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],tolA)
title('Average scaled data - tolA')
caxis manual
caxis ([bottom top])
colorbar
subplot(2,1,2)
imagesc(-2:2:10,[-1/2,1/2],f3)
title('Computational - tolA')
caxis manual
caxis ([bottom top])
colorbar

%%
%tolB kymographs

g3=spatialFRAP_tolB_d(d(1),d(2),d(3));
g3=g3./trapz(xs,g3(:,1));

bottom = min(min(min(tolB)),min(min(g3)));
top  = max(max(max(tolB)),max(max(g3)));

figure(6)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],tolB)
title('Average scaled data - tolB')
caxis manual
caxis ([bottom top])
colorbar
subplot(2,1,2)
imagesc(-2:2:10,[-1/2,1/2],g3)
title('Computational - tolB')
caxis manual
caxis ([bottom top])
colorbar

%%

%find average length of cells in um
lngth=cellfun('size',B_d.cells,1);
Ld=median(lngth)*B_d.pixelsize;
lngth=cellfun('size',B_nd.cells,1);
Lnd=median(lngth)*B_nd.pixelsize;

%find prebleach steady state solution
[w1,w2,w3,w4]=steady_state_d(d(1),d(2),d(3),Ld);
[w10,w20,w30,w40]=steady_state_nd(d(1),d(2),d(3),Lnd);

%%

%plot TolB concentration in inner and outer for dividing and non-dividing
%cells
x=-1/2:0.005:1/2;
figure(7)
clf
plot(x,w1,'color',[0.85, 0.325, 0.098],'DisplayName','Dividing - outer')
hold on
plot(x,w2,'color',[0, 0.447, 0.741],'DisplayName','Dividing - inner')
plot(x,w10,'--','color',[0.85, 0.325, 0.098],'DisplayName','Non-dividing - outer')
plot(x,w20,'--','color',[0, 0.447, 0.741],'Displayname','Non-dividing - inner')
hold off
legend
ylim([0 20000])
title('TolB concentration profile')

trapz(w1)
trapz(w2)

%%

%plot Pal concentration in inner and outer for dividing and non-dividing
%cells
x=-1/2:0.005:1/2;
figure(9)
clf
plot(x,w1,'color',[0.85, 0.325, 0.098],'DisplayName','Dividing - complex')
hold on
plot(x,w3,'color',[0, 0.447, 0.741],'DisplayName','Dividing - free')
plot(x,w4,'color',[0.9290, 0.6940, 0.1250],'DisplayName','Dividing - bound')
plot(x,w10,'--','color',[0.85, 0.325, 0.098],'DisplayName','Non-dividing - complex')
plot(x,w30,'--','color',[0, 0.447, 0.741],'Displayname','Non-dividing - free')
plot(x,w30,'--','color',[0.9290, 0.6940, 0.1250],'Displayname','Non-dividing - bound')
hold off
legend
%ylim([0 18000])
title('Pal concentration profile')

trapz(x,w1+w3+w4);
trapz(x,w10+w30+w40);

%%
%find Pal effective difusion coefficients

%find average length of cells in um
lngth=cellfun('size',P_d.cells,1);
Pd_lngth=median(lngth)/50;
lngth=cellfun('size',P_nd.cells,1);
Pnd_lngth=median(lngth)/50;
lngth=cellfun('size',A.cells,1);
PA_lngth=median(lngth)/50;
lngth=cellfun('size',B.cells,1);
PB_lngth=median(lngth)/50;

%fit to find effective diffusion constant from experimental data and
%adjust for cell length and pixelsize
Pnd_data=P_nd.avg;
[d2]=fitkymo(P_nd.t,Pnd_data,0.01);
Pnd_diff=d2*(Pnd_lngth*P_nd.pixelsize)^2;
Pd_data=P_d.avg;
[d2]=fitkymo(P_d.t,Pd_data,0.01);
Pd_diff=d2*(Pd_lngth*P_d.pixelsize)^2;
PA_data=A.avg;
[d2]=fitkymo(A.t,PA_data,0.01);
PA_diff=d2*(PA_lngth*A.pixelsize)^2;
PB_data=B.avg;
[d2]=fitkymo(B.t,PB_data,0.01);
PB_diff=d2*(PB_lngth*B.pixelsize)^2;

%adjust by initial frame shape to get deff
deff_d=Pd_diff./Pd_data(:,1)/length(Pd_data(:,1));
deff_nd=Pnd_diff./Pnd_data(:,1)/length(Pnd_data(:,1));
deff_A=PA_diff./PA_data(:,1)/length(PA_data(:,1));
deff_B=PB_diff./PB_data(:,1)/length(PB_data(:,1));
deff_d=interp1(-1/2:0.02:1/2,deff_d,-1/2:0.005:1/2)';
deff_nd=interp1(-1/2:0.02:1/2,deff_nd,-1/2:0.005:1/2)';
deff_A=interp1(-1/2:0.02:1/2,deff_A,-1/2:0.005:1/2)';
deff_B=interp1(-1/2:0.02:1/2,deff_B,-1/2:0.005:1/2)';

%find model deff
Deff=effective_diff_d(d(1),d(2),d(3));
Deff0=effective_diff_nd(d(1),d(2),d(3));
DeffA=effective_diff_A(d(1),d(2),0);%tolA mutant, no sink
DeffB=effective_diff_B(d(1),d(2),d(3));

%%
%plot model and experimental effective diffusion coefficients

x=-1/2:0.005:1/2;
figure(8)
clf
plot(x,deff_nd,'--','Color',[0, 0.4470, 0.7410],'DisplayName','Exp - non-dividing')
hold on
plot(x,deff_d,'--','Color',[0.8500, 0.3250, 0.0980],'DisplayName','Exp - dividing')
plot(x,deff_A,'--','Color',[0.9290, 0.6940, 0.1250],'DisplayName','Exp - tolA')
plot(x,deff_B,'--','Color',[0.4940, 0.1840, 0.5560],'DisplayName','Exp - tolB')
plot(x,Deff0,'Color',[0, 0.4470, 0.7410],'DisplayName','Non-dividing')
plot(x,Deff,'Color',[0.8500, 0.3250, 0.0980],'DisplayName','Dividing')
plot(x,DeffA,'Color',[0.9290, 0.6940, 0.1250],'DisplayName','tolA')
plot(x,DeffB,'Color',[0.4940, 0.1840, 0.5560],'DisplayName','tolB')
hold off
ylim([0 2e-3])
legend

Deff_10=Deff;

save('TolA_overexpression.mat','-append','Deff_10')
