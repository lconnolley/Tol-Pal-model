function [d,fval]=fitkymo4(tolb_d,pal_d,tolb_nd,pal_nd,guess)

tolb_d=interp1(-1/2:0.02:1/2,tolb_d,-1/2:0.01:1/2);
pal_d=interp1(-1/2:0.02:1/2,pal_d,-1/2:0.01:1/2);
tolb_nd=interp1(-1/2:0.02:1/2,tolb_nd,-1/2:0.01:1/2);
pal_nd=interp1(-1/2:0.02:1/2,pal_nd,-1/2:0.01:1/2);

lb=[0 0 1e7 1e4 0];
ub=[0.1 0.1 1e10 1e6 1e-2];

fun=@(d) cost(tolb_d(:,2:end), spatialFRAP_TolB_d(d(1),d(2),d(3),d(4),d(5)), pal_d(:,2:end), spatialFRAP_Pal_d(d(1),d(2),d(3),d(4),d(5)), tolb_nd(:,2:end), spatialFRAP_TolB_nd(d(1),d(2),d(3),d(4),d(5)), pal_nd(:,2:end), spatialFRAP_Pal_nd(d(1),d(2),d(3),d(4),d(5)));

%Pattern search
options=optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,'InitialMeshSize',0.1,'Display','iter');
[d,fval,exitflag,output] = patternsearch(fun,guess,[],[],[],[],lb,ub,[],options);

%Surrogate optimisation
%options = optimoptions('surrogateopt','UseParallel',true,'MaxFunctionEvaluations',1000);
%[d,fval,exitflag,output] = surrogateopt(fun,lb,ub,options);

%Multistart
%rng default %for reproducibility
%options = optimoptions(@fmincon,'UseParallel',true,'Display','iter');
%problem = createOptimProblem('fmincon','objective',fun,'lb',lb,'ub',ub,'x0',guess,'options',options);
%ms = MultiStart();
%[d,fval] = run(ms,problem,20);

%------------------------------------------------------------------

function out=cost(A,B,C,D,E,F,G,H)
Anorm = A / trapz(tolb_d(:,1));
Bnorm = B / trapz(B(:,1));
Cnorm = C / trapz(pal_d(:,1));
Dnorm = D / trapz(D(:,1));
Enorm = E / trapz(tolb_nd(:,1));
Fnorm = F / trapz(F(:,1));
Gnorm = G / trapz(pal_nd(:,1));
Hnorm = H / trapz(H(:,1));

out1=1e5*immsre(Anorm,Bnorm);
out2=1e5*immsre(Cnorm,Dnorm);
out3=1e5*immsre(Enorm,Fnorm);
out4=1e5*immsre(Gnorm,Hnorm);
out=out1+out2+out3+out4;

end

end