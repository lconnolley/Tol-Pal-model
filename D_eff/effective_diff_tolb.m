%effective diffusion coefficients for TolB from model

clear all

%%
%load data

Bd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_dividing_nopeaks.mat');
Bnd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_nondiv_2s.mat');

%%
%calculate average cell length in um

lngth=cellfun('size',Bd.cells,1);
Ld=median(lngth)*Bd.pixelsize;
lngth=cellfun('size',Bnd.cells,1);
Lnd=median(lngth)*Bnd.pixelsize;
%%
%set parameter choice

%From fitting, [Dc, Db, beta0, N]
d=[0.0129, 0.0057, 2.1748e9, 1.7e5];

Dc=d(1);
Db=d(2);
Df=Dc;              %Victor's paper
Dp=0.000;
alpha=5.4e4;        %Papadakos paper
beta0=d(3);
gamma=0.006;        %Papadakos paper
kon=1e-3;            %estimate was 1e5-1e6
koff=1;            %estimate was 1-10
N=d(4);
sigma=0.05;

%%
%generate model kymograph data

[w1,w2,w3,w4]=steady_state_d(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma,Ld);
[w10,w20,w30,w40]=steady_state_nd(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma,Lnd);

tolb=spatialFRAP_TolB_d(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma);
tolb0=spatialFRAP_TolB_nd(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma);

%%
%combine steady state and post bleach

data=[(w1+w2).*Bd.factor; tolb]';
data0=[(w10+w20).*Bnd.factor; tolb0]';

%%
%find deff

[Deff,q]=deffective(Bd.t,data);
[Deff0,d]=deffective(Bnd.t,data0);

%%
%calculate experimental deff

Bnd_data=Bnd.avg;
[d2]=fitkymo(Bnd.t,Bnd_data,0.01);
Bnd_diff=d2*(Lnd/50)^2;
Bd_data=Bd.avg;
[d2]=fitkymo(Bd.t,Bd_data,0.01);
Bd_diff=d2*(Ld/50)^2;

deff_d=Bd_diff./Bd_data(:,1)/length(Bd_data(:,1));
deff_nd=Bnd_diff./Bnd_data(:,1)/length(Bnd_data(:,1));
%%
%plot experimental and model deff

figure(3)
clf
plot(-1/2:0.02:1/2,deff_nd,'DisplayName','Exp - non-div')
hold on
plot(-1/2:0.02:1/2,deff_d,'DisplayName','Exp - div')
plot(-1/2:0.005:1/2,Deff0*(Lnd*0.005)^2,'DisplayName','Sim - non-div')
plot(-1/2:0.005:1/2,Deff*(Ld*0.005)^2,'DisplayName','Sim - div')
hold off
xlabel('Relative position')
ylabel('Effective diffusion coefficient')
%ylim([0 0.07])
legend

%%
%plot TolB concentration profile

%{
figure(14)
clf
plot(x,w1,'color',[0, 0.447, 0.741],'DisplayName','Dividing - outer')
hold on
plot(x,w2,'color',[0.85, 0.325, 0.098],'DisplayName','Dividing - inner')
plot(x,w10,'--','color',[0, 0.447, 0.741],'DisplayName','Non-dividing - outer')
plot(x,w20,'--','color',[0.85, 0.325, 0.098],'Displayname','Non-dividing - inner')
hold off
legend
title('TolB concentration profile')
%}
