function Deff=effective_diff_B(a,b,beta0)

B=load('../Import/tolB_dividing.mat');

pal=spatialFRAP_tolB_d(a,b,beta0);

guess=0.01;
d=fitkymo(B.t,pal,guess);
Deff=d./pal(:,1)/length(pal(:,1));

lngth=cellfun('size',B.cells,1);
L=median(lngth)*B.pixelsize;

Deff=Deff*(L*0.005)^2;

end