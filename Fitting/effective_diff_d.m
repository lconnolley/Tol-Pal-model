function Deff=effective_diff_d(a,b,beta0)

P=load('../Import/Pal_dividing_30s.mat');%Pal_dividing.mat');

pal=spatialFRAP_Pal_d(a,b,beta0);

guess=0.01;
d=fitkymo(P.t,pal,guess);
Deff=d./pal(:,1)/length(pal(:,1));

lngth=cellfun('size',P.cells,1);
L=median(lngth)*P.pixelsize;

Deff=Deff*(L*0.005)^2;

end