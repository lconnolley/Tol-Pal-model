clear all

%%
t0=0:60:60*60;
t=0:0.003:0.03;

Dc=0.1;
Db=0.5;
alpha=50;
beta0=100;

d=Dc/Db;
a=alpha/Db;
b=beta0/alpha;

d=50;
a=50;
b=1;

%%
[c0,b0]=pre_bleach(d,a,b,t0);

%%
[tran,tran2]=post_bleach(d,a,b,t);
[no_tran,no_tran2]=post_bleach(d,0,0,t);

fact=trapz(tran(1,:)+tran2(1,:))/trapz(c0(end,:)+b0(end,:));

T = [(c0(end,:)+b0(end,:)).*fact; tran+tran2]';
NT = [(c0(end,:)+b0(end,:)).*fact; no_tran+no_tran2]';

t0=diff(t);
t=[-t0(1); t']';
figure(1)
clf
subplot(2,1,1)
imagesc(t,[-1/2:1/2],T)
title('Simulated kymograph for homogenous transport')
subplot(2,1,2)
imagesc(t,[-1/2,1/2],NT)
title('Simulated kymograph for no transport')

trapz(T)+trapz(NT)


