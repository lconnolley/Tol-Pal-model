clear all;

B_d=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_dividing_5s_nopeaks.mat');
P_d=load('Pal_dividing.mat');

B_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_nondividing_5s.mat');
P_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/non-dividing.mat');

%%
%averaged cell

guess=[0.004 0.002 7.0e8 3e5 1.0e-4]; %a, b, beta0, N, kon

%renormalise data to tolb/pal concentration after bleach, change so 
%normalised to 1 on x, multiply by average concentration of total tolb/pal 
%in a cell (uM), multiply by area under bleached curve
tolb_d=B_d.avg*length(B_d.avg)*9.96e3*B_d.factor;
pal_d=P_d.avg*length(P_d.avg)*99.6e3*P_d.factor;
tolb_nd=B_nd.avg*length(B_nd.avg)*9.96e3*B_nd.factor;
pal_nd=P_nd.avg*length(P_nd.avg)*99.6e3*P_nd.factor;

%%
tic
[d,fval]=fitkymo4(tolb_d,pal_d,tolb_nd,pal_nd,guess);
toc
d          
fval
%%

%find prebleach steady state solution
[w1,w2,w3,w4]=steady_state_d(d(1),d(2),d(3),d(4),d(5));
[w10,w20,w30,w40]=steady_state_nd(d(1),d(2),d(3),d(4),d(5));

figure(3)
clf
subplot(2,1,1)
imagesc(-5:5:50,[-1/2,1/2],tolb_d)
title('Average scaled data - TolB (dividing)')
subplot(2,1,2)
f1=spatialFRAP_TolB_d(d(1),d(2),d(3),d(4),d(5));
imagesc(-5:5:50,[-1/2,1/2],[(w1'+w2')*B_d.factor,f1])
title('Computational - TolB (dividing)')

trapz(f1);



figure(4)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],pal_d)
title('Average scaled data - Pal (dividing)')
subplot(2,1,2)
g1=spatialFRAP_Pal_d(d(1),d(2),d(3),d(4),d(5));
imagesc(-2:2:10,[-1/2,1/2],[(w1'+w3'+w4')*P_d.factor,g1])
title('Computational - Pal (dividing)')

trapz([(w1'+w3'+w4')*P_d.factor,g1])

figure(5)
clf
subplot(2,1,1)
imagesc(-5:5:100,[-1/2,1/2],tolb_nd)
title('Average scaled data - TolB (non-dividing)')
subplot(2,1,2)
f2=spatialFRAP_TolB_nd(d(1),d(2),d(3),d(4),d(5));
imagesc(-5:5:100,[-1/2,1/2],[(w10'+w20')*B_nd.factor,f2])
title('Computational - TolB (non-dividing')

trapz(f2);

figure(6)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],pal_nd)
title('Average scaled data - Pal (non-dividing)')
subplot(2,1,2)
g2=spatialFRAP_Pal_nd(d(1),d(2),d(3),d(4),d(5));
imagesc(-2:2:10,[-1/2,1/2],[(w10'+w30'+w40')*P_nd.factor,g2])
title('Computational - Pal (non-dividing)')

trapz(g2);

