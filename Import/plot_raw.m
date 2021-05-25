 clear all
name='div2xbleach.mat';
load(name);

 t=[-30; t']';

i=1;
while i<length(cells)+1
    
    data=cells{i};
    %data=bindata(data,binfact);
    %data=data(3:end-2,:)*diag(1./sum(data(3:end-2,:),1));
    %npixels(i)=length(data(:,1));
    %data=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    
    figure(i)
    clf
    subplot(2,1,1)
    plot(data)
    xlim([1 51])
    
    subplot(2,1,2)
    %surf(t,-0.5:0.02:0.5,data,'EdgeColor','none');
    %view(2) 
    %xlim([-20, 300])
    imagesc(t,[-1/2,1/2],data)

    trapz(data)

    i=i+1;
end


