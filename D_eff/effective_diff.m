clear all

P_d=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Fitting/Pal_dividing.mat');
P_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/non-dividing.mat');

A=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolA_dividing.mat');
B=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolB_dividing.mat');

lngth=cellfun('length',P_d.cells);
Ld=mean(lngth)*P_d.pixelsize;
lngth=cellfun('length',P_nd.cells);
Lnd=mean(lngth)*P_nd.pixelsize;
lngth=cellfun('length',A.cells);
LA=mean(lngth)*A.pixelsize;
lngth=cellfun('length',B.cells);
LB=mean(lngth)*B.pixelsize;

%From fitting, [Dc, Db, beta0, N, kon
%d=[0.0102, 0.0056, 1.28e9, 1.62e5, 1.0e-2];
d=[0.0027, 0.0014, 7.95e8, 3e5, 1.0747e8];

Dc=d(1);            %Dc>0.0125
Db=d(2);            %Db<0.004
Df=Dc;              %Victor's paper
Dp=0.000;
alpha=5.4e4;        %Papadakos paper
beta0=d(3);
gamma=0.006;        %Papadakos paper
kon=d(5);            %estimate was 1e5-1e6
koff=10;            %estimate was 1-10
N=d(4);
sigma=0.05;

x=-1/2:0.005:1/2;

[w1,w2,w3,w4]=steady_state_d(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma,Ld);
[w10,w20,w30,w40]=steady_state_nd(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma,Lnd);
[w1A,w2A,w3A,w4A]=steady_state_d(Dc,Db,alpha,0,gamma,kon,koff,N,sigma,LA);
[w1B,w2B,w3B,w4B]=steady_state_d(Dc,Db,0,beta0,gamma,kon,koff,N,sigma,LB);

pal=spatialFRAP_Pal_d(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma);
pal0=spatialFRAP_Pal_nd(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma);
palA=spatialFRAP_Pal_A(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma);
palB=spatialFRAP_Pal_B(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma);

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



data=[(w1+w3+w4).*P_d.factor; pal]';
data0=[(w10+w30+w40).*P_nd.factor; pal0]';
dataB=[(w1B+w3B+w4B).*P_d.factor; palB]';
dataA=[(w1A+w3A+w4A).*P_d.factor; palA]';

Deff=deffective(P_d.t,data);
Deff0=deffective(P_nd.t,data0);
DeffB=deffective(B.t,dataB);
DeffA=deffective(A.t,dataA);

figure(11)
clf
plot(x,Deff0,'DisplayName','Non-dividing')
hold on
plot(x,Deff,'DisplayName','Dividing')
plot(x,DeffB,'DisplayName','tolB')
plot(x,DeffA,'DisplayName','tolA')
hold off
xlabel('Relative position')
ylabel('Effective diffusion coefficient')
legend


sol=spatialFRAP(P_d.t,pal(1,:),median(Deff/0.005^2));
bottom = min(min(min(pal)),min(min(sol)));
top  = max(max(max(pal)),max(max(sol)));

t1=-2:2:10;

figure(12)
clf
subplot(2,1,1)
imagesc(t1,[-1/2,1/2],data)
title('Averaged scaled data')
subplot(2,1,2)
imagesc(t1,[-1/2,1/2],[((w1+w3+w4).*P_d.factor)',sol])
title('Fokker-Planck')
caxis manual
caxis([bottom top]);

sol0=spatialFRAP(P_nd.t,pal0(1,:),median(Deff0/0.005^2));
bottom = min(min(min(pal0)),min(min(sol0)));
top  = max(max(max(pal0)),max(max(sol0)));

figure(13)
clf
subplot(2,1,1)
imagesc(t1,[-1/2,1/2],data0)
title('Average scaled data')
subplot(2,1,2)
imagesc(t1,[-1/2,1/2],[((w10+w20+w30).*P_nd.factor)',sol0])
title('Fokker-Planck')
caxis manual 
caxis([bottom top]);
