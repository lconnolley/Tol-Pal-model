clear

%%
%load data

z_1 = load('TolB_dividing_nopeaks.mat');
z_2 = load('TolB_nondiv_2s.mat');

%z_1 = load('Pal_dividing_30s.mat');
%z_2 = load('Pal_nondividing_30s.mat');

%%

conc_d=[];
lngth_d=[];

for i = 1:length(z_1.cells)
    
    data=z_1.cells{i};
    data=data(:,1);
    l=length(data)*z_1.pixelsize;
    conc_d=[conc_d; sum(data)/l];
    lngth_d=[lngth_d; l];
end

conc_nd=[];
lngth_nd=[];

for i = 1:length(z_2.cells)
    
    data=z_2.cells{i};
    data=data(:,1);
    l=length(data)*z_2.pixelsize;
    conc_nd=[conc_nd; sum(data)/l];
    lngth_nd=[lngth_nd; l];
end


figure(1)
clf
scatter(lngth_nd,conc_nd)
hold on
scatter(lngth_d,conc_d)
%ylim([0 1100])
xlabel('Cell length (um)')
ylabel('Total Fluorescence')

%%
%load data
clear all

z_1 = load('TolA_IPTG_distribution.mat');

%%

conc_d=[];
lngth_d=[];

for i = 1:84
    
    data=z_1.tolA_chr_d{i};
    l=length(data);
    conc_d=[conc_d; nansum(data)/l];
    lngth_d=[lngth_d; l*0.117];
end

conc_nd=[];
lngth_nd=[];

for i = 1:104
    
    data=z_1.tolA_chr_nd{i};
    l=length(data);
    conc_nd=[conc_nd; nansum(data)/l];
    lngth_nd=[lngth_nd; l*0.117];
end

figure(2)
clf
scatter(lngth_nd,conc_nd)
hold on
scatter(lngth_d,conc_d)
xlabel('Cell length (um)')
ylabel('Total Fluorescence')
ylim([0 750])

%correlation
[rho, pval] = corr([conc_nd; conc_d], [lngth_nd; lngth_d])



