%plot distribution of TolA in WT, KO, and arabinose induced
clear all

load('TolA_ara_distribution.mat')

%%
%sort data

tolA_KO=tolA_KO(~cellfun('isempty',tolA_KO));
tolA_chr=tolA_chr(~cellfun('isempty',tolA_chr));
tolA_2=tolA_2(~cellfun('isempty',tolA_2));
tolA_02=tolA_02(~cellfun('isempty',tolA_02));
tolA_002=tolA_002(~cellfun('isempty',tolA_002));
tolA_0002=tolA_0002(~cellfun('isempty',tolA_0002));
tolA_00002=tolA_00002(~cellfun('isempty',tolA_00002));
tolA_0=tolA_0(~cellfun('isempty',tolA_0));

%%
%find medians

avg_KO=[];
for i=1:length(tolA_KO)

    data=tolA_KO{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_KO=[avg_KO; a];

end
M_KO=median(avg_KO,1);

avg_chr=[];
for i=1:length(tolA_chr)

    data=tolA_chr{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_chr=[avg_chr; a];

end
M_chr=median(avg_chr,1);

avg_2=[];
for i=1:length(tolA_2)

    data=tolA_2{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_2=[avg_2; a];

end
M_2=median(avg_2,1);

avg_02=[];
for i=1:length(tolA_02)

    data=tolA_02{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_02=[avg_02; a];

end
M_02=median(avg_02,1);

avg_002=[];
for i=1:length(tolA_002)

    data=tolA_002{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_002=[avg_002; a];

end
M_002=median(avg_002,1);

avg_0002=[];
for i=1:length(tolA_0002)

    data=tolA_0002{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_0002=[avg_0002; a];

end
M_0002=median(avg_0002,1);

avg_00002=[];
for i=1:length(tolA_00002)

    data=tolA_00002{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_00002=[avg_00002; a];

end
M_00002=median(avg_00002,1);

avg_0=[];
for i=1:length(tolA_0)

    data=tolA_0{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_0=[avg_0; a];

end
M_0=median(avg_0,1);

%%
%plot figures

x=-1/2:0.02:1/2;

%KO and chromosomal TolA expression
figure(1)
clf
hold on
CI=bootci(1000,@nanmedian,avg_KO);
shadederror(x,nanmedian(avg_KO,1)',CI(1,:),CI(2,:),'KO')
CI=bootci(1000,@nanmedian,avg_chr);
shadederror(x,nanmedian(avg_chr,1)',CI(1,:),CI(2,:),'Chromsome')
hold off
xlabel('Relative position')
ylabel('Fluorescence')
legend

%KO and arabinose induced TolA expression
figure(2)
clf
hold on
CI=bootci(1000,@nanmedian,avg_chr);
shadederror(x,nanmedian(avg_chr,1)',CI(1,:),CI(2,:),'WT')
CI=bootci(1000,@nanmedian,avg_0);
shadederror(x,nanmedian(avg_0,1)',CI(1,:),CI(2,:),'0% ara')
CI=bootci(1000,@nanmedian,avg_002);
shadederror(x,nanmedian(avg_002,1)',CI(1,:),CI(2,:),'0.02% ara')
CI=bootci(1000,@nanmedian,avg_02);
shadederror(x,nanmedian(avg_02,1)',CI(1,:),CI(2,:),'0.2% ara')
hold off
xlabel('Relative position')
ylabel('Fluorescence')
legend

%plot median with background (KO) removed
M_chr=M_chr-M_KO;
M_chr(M_chr<0)=0;%set all negative values to 0
M_0=M_0-M_KO;
M_00002=M_00002-M_KO;
M_0002=M_0002-M_KO;
M_002=M_002-M_KO;
M_02=M_02-M_KO;
M_2=M_2-M_KO;

figure(3)
clf
plot(x,M_chr,'DisplayName','WT')
hold on
plot(x,M_0,'DisplayName','0% ara')
%plot(x,M_00002,'DisplayName','0.0002% ara')
%plot(x,M_0002,'DisplayName','0.002% ara')
plot(x,M_002,'DisplayName','0.02% ara')
plot(x,M_02,'DisplayName','0.2% ara')
%plot(x,M_2,'DisplayName','2% ara')
hold off
legend

%%
%finding sigma - width of TolA peak

%L=mean((cell_lengths(1,:)));
L=1;
x=-L/2:0.02*L:L/2;

figure(4)
clf
y=[zeros(18,1); M_chr(19:38)'; zeros(13,1)]';%take central peak only
i=trapz(x,y);
y=2/i * y;%normalise to 2

%shape of sink, beta
mu=0;
sigma=0.05*L;
beta=@(mu,x) normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma));%truncated normal
i = trapz(x,beta(mu,x));
beta=@(mu,x) 2/i*normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma));%normalise to 2

plot(x,y,'DisplayName','TolA')
hold on
plot(x,beta(mu,x),'DisplayName','Truncated normal')
xlim([min(x) max(x)])
legend

%%
%TolA overexpression background

z=M_02;
z(isnan(z))=760;
i=trapz(x,z);
z=2/i * z;%normalise to 2

%shape of sink, beta
mu=0;
background=20;
sigma=0.05*L;
beta=@(mu,x) normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma)) + background;%truncated normal
i = trapz(x,beta(mu,x));
beta=@(mu,x) 2/i*(normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma)) + background);%normalise to 2

figure(5)
clf
plot(x,z,'DisplayName','0.2% ara')
hold on
plot(x,beta(mu,x),'DisplayName','Normal + baseline')
hold off
legend


%%
%approx difference between peak and sides

WT=max(M_chr)/M_chr(:,12)
A0=max(M_0)/M_0(:,12)
A02=max(M_002)/M_002(:,12)
A2=max(M_02)/M_02(:,12)



