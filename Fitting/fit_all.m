clear all;

B_d=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_dividing_nopeaks.mat');
P_d=load('Pal_dividing.mat');

B_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_nondiv_2s_150s.mat');
P_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/non-dividing.mat');

A=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolA_dividing.mat');
B=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolB_dividing.mat');

%%
%averaged cell

guess=[0.0102 0.0056 1.28e9 1.62e5]; %Dc, Db, beta0, N

%renormalise data to tolb/pal concentration after bleach, change so 
%normalised to 1 on x, multiply by average concentration of total tolb/pal 
%in a cell (uM), multiply by area under bleached curve so same total
%fluorescene as after bleaching
tolb_d=B_d.avg*length(B_d.avg)*9.96e3*B_d.factor;
pal_d=P_d.avg*length(P_d.avg)*99.6e3*P_d.factor;
tolb_nd=B_nd.avg*length(B_nd.avg)*9.96e3*B_nd.factor;
pal_nd=P_nd.avg*length(P_nd.avg)*99.6e3*P_nd.factor;
tolA=A.avg*length(A.avg)*99.6e3*A.factor;   %beta0=0
tolB=B.avg*length(B.avg)*99.6e3*B.factor;   %alpha=0

tic
[d,fval]=fitkymo6(tolb_d,pal_d,tolb_nd,pal_nd,tolA,tolB,guess);
fval
toc

%find prebleach steady state solution
[w1,w2,w3,w4]=steady_state_d(d(1),d(2),d(3),d(4));
[w10,w20,w30,w40]=steady_state_nd(d(1),d(2),d(3),d(4));
[w1A,w2A,w3A,w4A]=steady_state_d(d(1),d(2),0,d(4));
[w1B,w2B,w3B,w4B]=steady_state_tolB(d(1),d(2),d(3),d(4));

figure(1)
clf
subplot(2,1,1)
imagesc(-5:5:50,[-1/2,1/2],tolb_d)
title('Average scaled data - TolB (dividing)')
subplot(2,1,2)
f1=spatialFRAP_TolB_d(d(1),d(2),d(3),d(4));
imagesc(-5:5:50,[-1/2,1/2],[(w1'+w3'+w4')*B_d.factor*0.1,f1])
title('Computational - TolB (dividing)')

trapz(f1)

figure(2)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],pal_d)
title('Average scaled data - Pal (dividing)')
subplot(2,1,2)
g2=spatialFRAP_Pal_d(d(1),d(2),d(3),d(4));
imagesc(-2:2:10,[-1/2,1/2],[(w1'+w2')*P_d.factor*10,g2])
title('Computational - Pal (dividing)')

trapz(g2)


figure(3)
clf
subplot(2,1,1)
imagesc(-5:5:100,[-1/2,1/2],tolb_nd)
title('Average scaled data - TolB (non-dividing)')
subplot(2,1,2)
f2=spatialFRAP_TolB_nd(d(1),d(2),d(3),d(4));
imagesc(-5:5:100,[-1/2,1/2],[(w10'+w30'+w40')*B_nd.factor*0.1,f2])
title('Computational - TolB (non-dividing)')

trapz(f2)

figure(4)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],pal_nd)
title('Average scaled data - Pal (non-dividing)')
subplot(2,1,2)
g2=spatialFRAP_Pal_nd(d(1),d(2),d(3),d(4));
imagesc(-2:2:10,[-1/2,1/2],[(w10'+w20')*P_nd.factor*10,g2])
title('Computational - Pal (non-dividing)')

trapz(g2)

figure(5)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],tolA)
title('Average scaled data - tolA')
subplot(2,1,2)
f3=spatialFRAP_tolA_d(d(1),d(2),0,d(4));
imagesc(-2:2:10,[-1/2,1/2],[(w1A'+w2A')*A.factor*10,f3])
title('Computational - tolA')

trapz(f3)

figure(6)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],tolB)
title('Average scaled data - tolB')
subplot(2,1,2)
g3=spatialFRAP_tolB_d(d(1),d(2),d(3),d(4));
imagesc(-2:2:10,[-1/2,1/2],[(w1B'+w2B')*B.factor*10,g3])
title('Computational - tolB')

trapz(g3)

