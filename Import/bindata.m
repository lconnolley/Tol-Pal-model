function out=bindata(y,n)
j=1;
for i=1:n:(size(y,1)-n+1)
    out(j,:)=mean(y(i:i+n-1,:),1);
    j=j+1;
end