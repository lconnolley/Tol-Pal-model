clear all;

%B_d=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_dividing_5s_no_peaks.mat');
P_d=load('Pal_dividing.mat');

%B_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_nondividing_5s.mat');
P_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/non-dividing.mat');
%%
%averaged cell

guess=[0.01 0.001 1.28e9 3e5]; %Dc, Db, beta0, N

%renormalise data to tolb/pal concentration after bleach, change so 
%normalised to 1 on x, multiply by average concentration of total tolb/pal 
%in a cell (uM), multiply by area under bleached curve
div=P_d.avg*length(P_d.avg)*99.6e3*P_d.factor;
non_div=P_nd.avg*length(P_nd.avg)*99.6e3*P_nd.factor;
tic
[d,fval]=fitkymo3(div,non_div,guess);
d           %*(binfact*pixelsize)^2
%fval
toc
%find prebleach steady state solution
[w1,w2,w3,w4]=steady_state_d(d(1),d(2),d(3),d(4));
[w10,w20,w30,w40]=steady_state_nd(d(1),d(2),d(3),d(4));

%%

figure(1)
clf
subplot(2,1,1) 
imagesc(-2:2:10,[-1/2,1/2],div)
title('Average scaled data - Pal dividing')
subplot(2,1,2)
f=spatialFRAP_Pal_d(d(1),d(2),d(3),d(4));
imagesc(-2:2:10,[-1/2,1/2],[(w1'+w3'+w4')*P_d.factor,f])
title('Computational - Pal dividing')

trapz(f);

figure(2)
clf
subplot(2,1,1)
imagesc(-2:2:10,[-1/2,1/2],non_div)
title('Average scaled data - Pal non-dividing')
subplot(2,1,2)
g=spatialFRAP_Pal_nd(d(1),d(2),d(3),d(4));
imagesc(-2:2:10,[-1/2,1/2],[(w10'+w30'+w40')*P_nd.factor,g])
title('Computational - Pal non-dividing')

trapz(g);


%plot TolB concentration profiles
x=-1/2:0.01:1/2;

figure(3)
clf
plot(x,w1,'color',[0, 0.447, 0.741],'DisplayName','Dividing - outer')
hold on
plot(x,w2,'color',[0.85, 0.325, 0.098],'DisplayName','Dividing - inner')
plot(x,w10,'--','color',[0, 0.447, 0.741],'DisplayName','Non-dividing - outer')
plot(x,w20,'--','color',[0.85, 0.325, 0.098],'Displayname','Non-dividing - inner')
hold off
legend
title('TolB concentration profile')


