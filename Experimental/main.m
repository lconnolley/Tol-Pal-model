%import excel data, process, and plot
clear all

%% give details of excel file

name='../Import/frap pal-mch WT dividing.xlsx';%excel filename
last=23;%number of pages in excel file
pxlsz=0.0976;%pixelsize of microscope used

t=0:30:600;%time fluorescence is recorded for
binfact=1;%adjust for oversampling from confocal microscope

mat_name='../Import/Pal_dividing_30s.mat';%name of .mat file to save to

%% if .mat file already exists, load this and plot data
if exist(mat_name,'file')
    disp('Using existing .mat file with this name')
    z=load(mat_name);
    [d]=plot_data(z.t,z.D_median,z.D,z.binfact,z.pixelsize,z.avg);

%if .mat file doesn't exist proceed with import and fitting
else
%% import data from excel

[pixelsize,cells]=importdata(name,last,pxlsz);

%% process data

[avg,D_median,D,results,bleach,factor]=fit_all_cells(t,pixelsize,binfact,cells);

%% plot data

[d]=plot_data(t,D_median,D,binfact,pixelsize,avg);

%% save to .mat file

save(mat_name,'t','pixelsize','cells','binfact','avg','D_median','D','results','bleach','factor')

end