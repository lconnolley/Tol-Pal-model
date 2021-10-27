function Deff=effective_diff_A(a,b,beta0)

A=load('../Import/tolA_dividing.mat');

pal=spatialFRAP_tolA_d(a,b,beta0);

guess=0.01;
d=fitkymo(A.t,pal,guess);
Deff=d./pal(:,1)/length(pal(:,1));

lngth=cellfun('size',A.cells,1);
L=median(lngth)*A.pixelsize;

Deff=Deff*(L*0.005)^2;

end