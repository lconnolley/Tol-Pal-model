clear all

for i=1:8
    
    data=xlsread('../Experimental/tolA distribution raw.xlsx',i);
    data=data(3:end,2:end);
    
    for j=1:size(data,2)%cell
        L=[];
        I=~isnan(data(:,j));
        L=find(I,1,'last');
        cells{i,j}(1:L,1)=data(1:L,j);
        cell_lengths(i,j)=(L*0.117)-0.117;
    end
end

tolA_chr_d=cells(1,:);
tolA_chr_nd=cells(2,:);
tolA_0_d=cells(3,:);
tolA_0_nd=cells(4,:);
tolA_05_d=cells(5,:);
tolA_05_nd=cells(6,:);
tolA_50_d=cells(7,:);
tolA_50_nd=cells(8,:);

save('../Import/TolA_ara_distribution.mat','tolA*','cell_lengths','cells')

