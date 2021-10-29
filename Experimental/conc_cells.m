% plot TolA and TolB correlation with length, figure 2 - S1

clear all

%% load data

z_1 = load('../Experimental/TolB_dividing.mat');
z_2 = load('../Experimental/TolB_nondiv_2s.mat');

%z_1 = load('../Import/Pal_dividing_30s.mat');
%z_2 = load('../Import/Pal_nondividing_30s.mat');

%%

conc_d=[];
lngth_d=[];

for i = 1:length(z_1.cells)
    
    data=z_1.cells{i};
    data=data(:,1);
    l=length(data)*z_1.pixelsize;
    conc_d=[conc_d; nansum(data)/l];
    lngth_d=[lngth_d; l];
end

conc_nd=[];
lngth_nd=[];

for i = 1:length(z_2.cells)
    
    data=z_2.cells{i};
    data=data(:,1);
    l=length(data)*z_2.pixelsize;
    conc_nd=[conc_nd; nansum(data)/l];
    lngth_nd=[lngth_nd; l];
end


figure(1)
clf
scatter(lngth_nd,conc_nd)
hold on
scatter(lngth_d,conc_d)
ylim([0 1100])
%ylim([0 3.5e5])
xlabel('Cell length (um)')
ylabel('Mean Fluorescence')

%correlation
[rho1, pval1] = corr([conc_nd; conc_d], [lngth_nd; lngth_d])


