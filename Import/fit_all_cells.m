function [avg,D_median,D,results,bleach,factor] = fit_all_cells(t,pixelsize,binfact,cells)
%process data imported via importdata.m

guess=1e-2/(binfact*pixelsize)^2; %in units of binned pixels
avg=0;
avg2=0;

for i=1:length(cells)
data=cells{i};
data=bindata(data,binfact);%use pixel binning to correct for oversampling (confocal microscope)

data2=data(4:end-3,:);%take off end, no normalisation(3:end-2)
data=data(4:end-3,:)*diag(1./sum(data(4:end-3,:),1));%take off the ends and normalise

npixels(i)=length(data(:,1));

avg=avg+interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
avg2=avg2+interp1(([1:npixels(i)]-1)/(npixels(i)-1),data2,0:0.02:1);

[d1,results(i,2),d2,results(i,4)]=fitkymo(t,data,guess); %fit
results(i,1)=d1*(binfact*pixelsize)^2;
results(i,3)=d2*(binfact*pixelsize)^2;

D_median(i)=median(results(i,3)./data(:,1)/length(data(:,1)));%%median (Fokker Planck) Deff accross each cell

D(:,i)=interp1(([1:npixels(i)]-1)/(npixels(i)-1),results(i,3)./data(:,1)/length(data(:,1)),0:0.02:1); %scaled Deff for all cells

end

avg=avg/length(cells);%average of all scaled cells
avg=avg*diag(1./sum(avg(1:end,:),1));

bleach=avg2(:,2)./avg2(:,1);%shape of bleached area
factor=trapz(0:0.02:1,bleach);

end