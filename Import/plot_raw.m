%plot raw kymographs and line profiles of every cell in file

clear

load('Pal_dividing_30s.mat');

 t=[-2; t']';

i=1;
while i<length(cells)+1
    
    data=cells{i};
    %remove edges and interpolate data 
    %data=bindata(data,binfact);
    %data=data(3:end-2,:)*diag(1./sum(data(3:end-2,:),1));
    %npixels(i)=length(data(:,1));
    %data=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    
    figure(i)
    clf
    subplot(2,1,1)
    plot(data)
    
    subplot(2,1,2)
    imagesc(t,[-1/2,1/2],data)

    trapz(data);

    i=i+1;
end


