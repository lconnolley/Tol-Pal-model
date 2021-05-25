clear all

z_d=load('Pal_dividing.mat');
z_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Models/non-dividing.mat');

%%
avg_nd=0;
A_nd=[];
for i=1:length(z_nd.cells)
data=z_nd.cells{i};

data=data(3:end-2,:)*diag(1./sum(data(3:end-2,:),1));%take off the ends and normalise

npixels(i)=length(data(:,1));
avg_nd=avg_nd+interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);

a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
a=a(:,1);
a=a';
a=a./trapz(a);
A_nd=[A_nd;a];

end

avg_d=0;
A_d=[];
for i=1:length(z_d.cells)
data=z_d.cells{i};

data=data(3:end-2,:)*diag(1./sum(data(3:end-2,:),1));%take off the ends and normalise

npixels(i)=length(data(:,1));
avg_d=avg_d+interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);

a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
a=a(:,1);
a=a';
a=a./trapz(a);
A_d=[A_d;a];

end

avg_nd=avg_nd./trapz(avg_nd);
avg_d=avg_d./trapz(avg_d);

figure(1)
clf
A_nd=A_nd';
CI=bootci(1000,@median,A_nd');
shadederror(0:0.02:1,median(A_nd,2)',CI(1,:),CI(2,:),'Non-Dividing')
hold on

A_d=A_d';
CI=bootci(1000,@median,A_d');
shadederror(0:0.02:1,median(A_d,2)',CI(1,:),CI(2,:),'Dividing')
hold off

legend
xlabel('Relative Position')
ylabel('Fluorescence')
ylim([0 0.04])



