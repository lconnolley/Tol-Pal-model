%plot experimental mutant kymographs
clear all

%%
%load data

z_A=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolA_dividing.mat');
z_AKO=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolAKO_dividing.mat');
z_B=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolB_dividing.mat');

%%
%process data into form needed

%renormalise data to tolb/pal concentration after bleach, change so 
%normalised to 1 on x, multiply by average concentration of total tolb/pal 
%in a cell (uM), multiply by area under bleached curve
A=z_A.avg*length(z_A.avg)*99.6e3;%*B_d.factor;
AKO=z_AKO.avg*length(z_AKO.avg)*99.6e3;
B=z_B.avg*length(z_B.avg)*99.6e3;

%%
%plot kymographs

figure(1)
clf
imagesc(-2:2:10,[-1/2,1/2],A)
title('tolA')

figure(2)
clf
imagesc(-2:2:10,[-1/2,1/2],AKO)
title('tolAKO')

figure(3)
clf
imagesc(-2:2:10,[-1/2,1/2],B)
title('tolB')

