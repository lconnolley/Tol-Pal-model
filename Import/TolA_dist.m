clear all

load('TolA_ara_distribution.mat')
%%
%{
%Trying to do everything within a for loop
avg=[];%create a ?x? matrix of NaN

for i=1:length(cells)
    tolA=cells(i,:);
    for j=1:length(tolA)
        data=tolA{i};
        npixels(i)=length(data(:,1));
        a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
        avg
%}

tolA_KO=tolA_KO(~cellfun('isempty',tolA_KO));
tolA_chr=tolA_chr(~cellfun('isempty',tolA_chr));
tolA_2=tolA_2(~cellfun('isempty',tolA_2));
tolA_02=tolA_02(~cellfun('isempty',tolA_02));
tolA_002=tolA_002(~cellfun('isempty',tolA_002));
tolA_0002=tolA_0002(~cellfun('isempty',tolA_0002));
tolA_00002=tolA_00002(~cellfun('isempty',tolA_00002));
tolA_0=tolA_0(~cellfun('isempty',tolA_0));

%--------------------------------------------------------------
avg_KO=[];
for i=1:length(tolA_KO)

    data=tolA_KO{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_KO=[avg_KO; a];

end
M_KO=mean(avg_KO,1);

%--------------------------------------------------------------
avg_chr=[];
for i=1:length(tolA_chr)

    data=tolA_chr{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_chr=[avg_chr; a];

end
M_chr=mean(avg_chr,1);

%--------------------------------------------------------------
avg_2=[];
for i=1:length(tolA_2)

    data=tolA_2{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_2=[avg_2; a];

end
M_2=mean(avg_2,1);

%--------------------------------------------------------------
avg_02=[];
for i=1:length(tolA_02)

    data=tolA_02{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_02=[avg_02; a];

end
M_02=mean(avg_02,1);

%--------------------------------------------------------------
avg_002=[];
for i=1:length(tolA_002)

    data=tolA_002{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_002=[avg_002; a];

end
M_002=mean(avg_002,1);

%--------------------------------------------------------------
avg_0002=[];
for i=1:length(tolA_0002)

    data=tolA_0002{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_0002=[avg_0002; a];

end
M_0002=mean(avg_0002,1);

%--------------------------------------------------------------
avg_00002=[];
for i=1:length(tolA_00002)

    data=tolA_00002{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_00002=[avg_00002; a];

end
M_00002=mean(avg_00002,1);

%--------------------------------------------------------------
avg_0=[];
for i=1:length(tolA_0)

    data=tolA_0{i};
    npixels(i)=length(data(:,1));
    a=interp1(([1:npixels(i)]-1)/(npixels(i)-1),data,0:0.02:1);
    avg_0=[avg_0; a];

end
M_0=mean(avg_0,1);


%figures-------------------------------------------------------

x=-1/2:0.02:1/2;

figure(1)
clf
hold on
plot(x,M_KO,'DisplayName','KO')
plot(x,M_chr,'DisplayName','Chromosome')
hold off
xlabel('Relative position')
ylabel('Fluorescence')
legend

figure(2)
clf
plot(x,M_KO,'DisplayName','KO')
hold on
plot(x,M_0,'DisplayName','0% ara')
plot(x,M_00002,'DisplayName','0.0002% ara')
plot(x,M_0002,'DisplayName','0.002% ara')
plot(x,M_002,'DisplayName','0.02% ara')
%plot(x,M_02,'DisplayName','0.2% ara')
%plot(x,M_2,'DisplayName','2% ara')
hold off
xlabel('Relative position')
ylabel('Fluorescence')
legend




