%import excel data, process, and plot

clear

%%
%give details of excel file

name='frap pal-mch pgfp-tolA 0.2 ara.xlsx';%excel filename
last=23;%number of pages in excel file
pxlsz=0.0976;%pixelsize of microscope used

t=0:30:600;%time fluorescence is recorded for
binfact=1;%adjust for oversampling from confocal microscope

mat_name='Pal_02_ara_remove_bright_cells.mat';%name of .mat file saved to

%%
%import data from excel
[pixelsize,cells]=importdata(name,last,pxlsz);

%%
%process data
[avg,D_median,D,results,bleach,factor]=fit_all_cells(t,pixelsize,binfact,cells);

%%
%plot data
[d]=plot_data(t,D_median,D,binfact,pixelsize,avg);

%%
%save to .mat file
save(mat_name,'t','pixelsize','cells','binfact','avg','D_median','D','results','bleach','factor')

