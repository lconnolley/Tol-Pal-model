function [pixelsize,cells]=importdata(name,last,pxlsz)

%import SpatialFRAP data from excel files

data2=xlsread(name,1);
background=data2(:,1:end);

for i=2:last %frame, adjust for respective data
        
data=xlsread(name,i);

data=data(1:end,1:end);

for j=1:size(data,2)%cell
    L=[];
    I=~isnan(data(:,j));
    L=find(I,1,'last');
    cells{j}(1:L,i-1)=data(1:L,j);
    cell_lengths(i-1,j)=L;
    
end
end

%find the median pixelsize used
pixelsize=median(((cell_lengths(2,:)*pxlsz)-pxlsz)./cell_lengths(2,:)');
pixelsize=median(pixelsize);

%remove background
for j=1:size(data,2)
    
     total(j)=mean(cells{j}(:,1))-background(1,j);
     cells{j}(cells{j}(:)<0)=0;
    
     %cells{j}=cells{j}*diag(1./sum(cells{j},1));%normalise cells
    
end

%remove any 'bad' cells
cells{:,17}=[];%remove cell 4
cells{:,19}=[];%remove cell 4
cells{:,21}=[];%remove cell 4
cells{:,22}=[];%remove cell 4
cells{:,27}=[];%remove cell 4
cells=cells(~cellfun('isempty',cells));

end
