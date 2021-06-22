clear all

for i=1:8
    
    data=xlsread('tolA ara distribution raw.xlsx',i);
    data=data(1:end,1:end);
    
    for j=1:size(data,2)%cell
        L=[];
        I=~isnan(data(:,j));
        L=find(I,1,'last');
        cells{i,j}(1:L,1)=data(1:L,j);
        cell_lengths(i,j)=(L*0.117)-0.117;
    end
end

tolA_chr=cells(1,:);
tolA_2=cells(2,:);
tolA_02=cells(3,:);
tolA_002=cells(4,:);
tolA_0002=cells(5,:);
tolA_00002=cells(6,:);
tolA_0=cells(7,:);
tolA_KO=cells(8,:);

%save('TolA_ara_distribution.mat','tolA*','cell_lengths','cells')

