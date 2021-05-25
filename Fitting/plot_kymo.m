clear all;

B_d=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_dividing_nopeaks.mat');
P_d=load('Pal_dividing.mat');

B_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_nondiv_2s_150s.mat');
P_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/non-dividing.mat');

A=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolA_dividing.mat');
B=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolB_dividing.mat');

%%
%averaged cell

%d from fitting, [Dc, Db, beta0, N];
%d=[0.0999 0.001 7.0e8 3.0e5];%from fitting
d=[0.0102 0.0056 1.28e9 1.62e5];

%renormalise data to tolb/pal concentration after bleach, change so 
%normalised to 1 on x, multiply by average concentration of total tolb/pal 
%in a cell (uM), multiply by area under bleached curve
tolb_d=B_d.avg*length(B_d.avg)*9.96e3*B_d.factor;
pal_d=P_d.avg*length(P_d.avg)*99.6e3*P_d.factor;
tolb_nd=B_nd.avg*length(B_nd.avg)*9.96e3*B_nd.factor;
pal_nd=P_nd.avg*length(P_nd.avg)*99.6e3*P_nd.factor;
tolA=A.avg*length(A.avg)*99.6e3*A.factor;   %beta0=0
tolB=B.avg*length(B.avg)*99.6e3*B.factor;   %alpha=0

%find prebleach steady state solution
[w1,w2,w3,w4]=steady_state_d(d(1),d(2),d(3),d(4));
[w10,w20,w30,w40]=steady_state_nd(d(1),d(2),d(3),d(4));
[w1A,w2A,w3A,w4A]=steady_state_d(d(1),d(2),0,d(4));
[w1B,w2B,w3B,w4B]=steady_state_tolB(d(1),d(2),d(3),d(4));


figure(11)
clf
subplot(2,1,1)
imagesc(-2:2:150,[-1/2,1/2],tolb_d)
title('Average scaled data - TolB (dividing)')
subplot(2,1,2)
f1=spatialFRAP_TolB_d(d(1),d(2),d(3),d(4));
imagesc(-2:2:150,[-1/2,1/2],[(w1'+w2')*B_d.factor*0.9170,f1])
title('Computational - TolB (dividing)')

TBd=trapz(tolb_d);

trapz([(w1'+w2')*B_d.factor*0.9170,f1]);

figure(12)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],pal_d)
title('Average scaled data - Pal (dividing)')
subplot(2,1,2)
g1=spatialFRAP_Pal_d(d(1),d(2),d(3),d(4));
imagesc(-2:2:10,[-1/2,1/2],[(w1'+w3'+w4')*P_d.factor,g1])
title('Computational - Pal (dividing)')

Pd=trapz(pal_d);

trapz([(w1'+w3'+w4')*P_d.factor,g1]);

figure(13)
clf
subplot(2,1,1)
imagesc(-2:2:150,[-1/2,1/2],tolb_nd)
title('Average scaled data - TolB (non-dividing)')
subplot(2,1,2)
f2=spatialFRAP_TolB_nd(d(1),d(2),d(3),d(4));
imagesc(-2:2:150,[-1/2,1/2],[(w10'+w20')*B_nd.factor,f2])
title('Computational - TolB (non-dividing)')

TB_nd=trapz(tolb_nd);

trapz([(w10'+w20')*B_nd.factor,f2]);

figure(14)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],pal_nd)
title('Average scaled data - Pal (non-dividing)')
subplot(2,1,2)
g2=spatialFRAP_Pal_nd(d(1),d(2),d(3),d(4));
imagesc(-2:2:10,[-1/2,1/2],[(w10'+w30'+w40')*P_nd.factor,g2])
title('Computational - Pal (non-dividing)')

Pnd=trapz(pal_nd);
trapz([(w10'+w30'+w40')*P_nd.factor,g2]);

figure(15)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],tolA)
title('Average scaled data - tolA')
subplot(2,1,2)
f3=spatialFRAP_tolA_d(d(1),d(2),0,d(4));
imagesc(-2:2:10,[-1/2,1/2],[(w1A'+w3A'+w4A')*A.factor,f3])
title('Computational - tolA')

trapz([(w1A'+w3A'+w4A')*A.factor,f3]);

figure(16)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],tolB)
title('Average scaled data - tolB')
subplot(2,1,2)
g3=spatialFRAP_tolB_d(d(1),d(2),d(3),d(4));
imagesc(-2:2:10,[-1/2,1/2],[(w1B'+w3B'+w4B')*B.factor,g3])
title('Computational - tolB')

trapz([(w1B'+w3B'+w4B')*B.factor,g3]);

x=-1/2:0.01:1/2;
figure(17)
clf
plot(x,w1,'color',[0, 0.447, 0.741],'DisplayName','Dividing - outer')
hold on
plot(x,w2,'color',[0.85, 0.325, 0.098],'DisplayName','Dividing - inner')
plot(x,w10,'--','color',[0, 0.447, 0.741],'DisplayName','Non-dividing - outer')
plot(x,w20,'--','color',[0.85, 0.325, 0.098],'Displayname','Non-dividing - inner')
hold off
legend
title('TolB concentration profile')

