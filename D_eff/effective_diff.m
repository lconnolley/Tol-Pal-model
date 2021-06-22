%effective diffusion coefficients for Pal from model

clear all

%%
%load data

P_d=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Fitting/Pal_dividing.mat');
P_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/non-dividing.mat');

A=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolA_dividing.mat');
B=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolB_dividing.mat');
%%
%calculate average cell length in um

lngth=cellfun('size',P_d.cells,1);
Ld=median(lngth)*P_d.pixelsize;
lngth=cellfun('size',P_nd.cells,1);
Lnd=median(lngth)*P_nd.pixelsize;
lngth=cellfun('size',A.cells,1);
LA=median(lngth)*A.pixelsize;
lngth=cellfun('size',B.cells,1);
LB=median(lngth)*B.pixelsize;
%%
%set parameter choice

%From fitting, [Dc, Db, beta0, N]
d=[0.0289, 0.0067, 5.6685e10, 1.7e5];

Dc=d(1);
Db=d(2);
Df=Dc;                  %Victor's paper
Dp=0.000;
alpha=5.4e4;            %Papadakos paper
beta0=d(3);
gamma=0.006;            %Papadakos paper
kon=1e-3;               %estimate was 1e5-1e6
koff=1;                 %estimate was 1-10
N=d(4);
sigma=0.05;

%%
%generate model kymograph data

[w1,w2,w3,w4]=steady_state_d(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma,Ld);
[w10,w20,w30,w40]=steady_state_nd(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma,Lnd);
[w1A,w2A,w3A,w4A]=steady_state_d(Dc,Db,alpha,0,gamma,kon,koff,N,sigma,LA);
[w1B,w2B,w3B,w4B]=steady_state_d(Dc,Db,0,beta0,gamma,kon,koff,N,sigma,LB);

pal=spatialFRAP_Pal_d(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma);
pal0=spatialFRAP_Pal_nd(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma);
palA=spatialFRAP_Pal_A(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma);
palB=spatialFRAP_Pal_B(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma);

%%
%plot concentration profiles of Pal
%{
figure(10)
clf
subplot(3,1,1)
plot(x,w10+w30+w40,'DisplayName','Non-dividing')
hold on
plot(x,w1+w3+w4,'DisplayName','Dividing')
hold off
title('Pal pre bleach')
subplot(3,1,2)
plot(x,pal0(1,:),'DisplayName','Non-dividing')
hold on
plot(x,pal(1,:),'DisplayName','Dividing')
hold off
title('Pal after bleaching')
legend
subplot(3,1,3)
plot(x,pal0(end,:),'DisplayName','Non-dividing')
hold on
plot(x,pal(end,:),'DisplayName','Dividing')
hold off
title('Pal after recovery')
legend
%}

%plot concentration profile of TolB
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

%%
%combine steady state and post bleach

data=[(w1+w3+w4).*P_d.factor; pal]';
data0=[(w10+w30+w40).*P_nd.factor; pal0]';
dataB=[(w1B+w3B+w4B).*P_d.factor; palB]';
dataA=[(w1A+w3A+w4A).*P_d.factor; palA]';

%%

%find deff
Deff=deffective(P_d.t,data);
Deff0=deffective(P_nd.t,data0);
DeffB=deffective(B.t,dataB);
DeffA=deffective(A.t,dataA);

%%

x=-1/2:0.005:1/2;

%plot effective diffusion coefficient
figure(1)
clf
plot(x,Deff0*(Lnd*0.005)^2,'DisplayName','Non-dividing','LineWidth',2)
hold on
plot(x,Deff*(Ld*0.005)^2,'DisplayName','Dividing','LineWidth',2)
plot(x,DeffB*(LB*0.005)^2,'DisplayName','tolB','LineWidth',2)
plot(x,DeffA*(LA*0.005)^2,'DisplayName','tolA','LineWidth',2)
hold off
xlabel('Relative position')
ylabel('Effective diffusion coefficient')
%ylim([0 1e-3])
legend
