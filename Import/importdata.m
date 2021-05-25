clear all;

t=0:30:600;%adjust for respective data
binfact=1;

data2=xlsread('nondividing 2xbleach.xlsx',1);
background=data2(:,1:end);

for i=2:23 %frame, adjust for respective data
        
data=xlsread('nondividing 2xbleach.xlsx',i);

data=data(1:end,1:end);

for j=1:size(data,2)%cell
    L=[];
    I=~isnan(data(:,j));
    L=find(I,1,'last');
    cells{j}(1:L,i-1)=data(1:L,j);
    cell_lengths(i-1,j)=L;
    
end
end


pixelsize=median(((cell_lengths(2,:)*0.0976)-0.0976)./cell_lengths(2,:)');
pixelsize=median(pixelsize);

for j=1:size(data,2)
    
     total(j)=mean(cells{j}(:,1))-background(1,j);
     cells{j}(cells{j}(:)<0)=0;
    
     %cells{j}=cells{j}*diag(1./sum(cells{j},1));%normalise cells
    
end

%cells{:,4}=[];%remove cell 4
%cells=cells(~cellfun('isempty',cells));

save('nondiv2xbleach.mat','t','pixelsize','cells','binfact')


%{
% data2=xlsread('background+pole dividing cells distribution raw.xlsx',1);
% data2=data2(2:end,4);
% data2=data2(~isnan(data2));
% for i=1:30
%     x=[];
%     %cols=1+[i,i+30,i+2*30,i+3*30,i+4*30,i+5*30,i+6*30];
%     bk=data2((i-1)*13+[1:13]);
%     bk=bk([1:2:13]);%every second frame
%     background(:,i)=bk;
% end
% 
% 
% for i=3:9
% data=xlsread('background+pole dividing cells distribution raw.xlsx',i);
% for j=1:size(data,2)
%     L=[];
%     I=~isnan(data(:,j));
%     L=find(I,1,'last');
%     cells{j}(1:L,i-2)=data(1:L,j);
%     cell_lengths(i-2,j)=L;
% end
% end
% 
% for j=1:size(data,2)
%     
%      total(j)=mean(cells{j}(:,1))-background(1,j);
%      cells{j}(cells{j}(:)<0)=0;
% 
%     cells{j}=cells{j}*diag(1./sum(cells{j},1));
%     
%     
% end
% 
% save('dividing_pole.mat','cells')
%}