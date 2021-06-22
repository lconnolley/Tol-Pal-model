function [d,fval]=fitkymo_pal(div,nondiv,deff_d,deff_nd,guess)

div=interp1(-1/2:0.02:1/2,div,-1/2:0.005:1/2);
nondiv=interp1(-1/2:0.02:1/2,nondiv,-1/2:0.005:1/2);
deff_d=interp1(-1/2:0.02:1/2,deff_d,-1/2:0.005:1/2)';
deff_nd=interp1(-1/2:0.02:1/2,deff_nd,-1/2:0.005:1/2)';

fun=@(d) cost(div, spatialFRAP_Pal_d(d(1),d(2),d(3)), nondiv, spatialFRAP_Pal_nd(d(1),d(2),d(3)), deff_d, effective_diff_d(d(1),d(2),d(3)), deff_nd, effective_diff_nd(d(1),d(2),d(3)));
%Pattern search
options=optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,'InitialMeshSize',0.1,'Display','iter');%,'StepTolerance',1e-9);
[d,fval,~,~] = patternsearch(fun,guess,[],[],[],[],[0 1 1e7],[1 10 1e11],[],options); %a>0 s.t. Dc>Db, b>1 s.t. Db,Dc>0

%------------------------------------------------------------------

function out=cost(A,B,C,D,E,F,G,H)
Anorm = A / trapz(div(:,1));
Bnorm = B / trapz(B(:,1));
Cnorm = C / trapz(nondiv(:,1));
Dnorm = D / trapz(D(:,1));

out1=1e5*immse(Anorm,Bnorm);
out2=1e5*immse(Cnorm,Dnorm);

out3=1e8*immse(E,F);
out4=1e8*immse(G,H);

out=out1+out2+out3+out4;

end

end