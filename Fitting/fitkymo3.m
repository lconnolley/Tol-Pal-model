function [d,fval]=fitkymo3(div,nondiv,guess)

div=interp1(-1/2:0.02:1/2,div,-1/2:0.01:1/2);
nondiv=interp1(-1/2:0.02:1/2,nondiv,-1/2:0.01:1/2);

fun=@(d) cost(div(:,2:end), spatialFRAP_Pal_d(d(1),d(2),d(3),d(4)), nondiv(:,2:end), spatialFRAP_Pal_nd(d(1),d(2),d(3),d(4)));
%Pattern search
options=optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,'InitialMeshSize',0.1,'Display','iter');%,'StepTolerance',1e-9);
[d,fval,exitflag,output] = patternsearch(fun,guess,[],[],[],[],[0 0 1e7 1e4],[0.05 0.01 1e10 1e6],[],options); %Dc, Db, beta0, N, 

%Surrogate optimisation
%options = optimoptions('surrogateopt','UseParallel',true,'MaxFunctionEvaluations',1000);
%[d,fval,exitflag,output] = surrogateopt(fun,[0 0 0],[10 1e20 1e20],options);

%Multistart
%rng default %for reproducibility
%options = optimoptions(@fmincon,'UseParallel',true,'Display','iter');
%problem = createOptimProblem('fmincon','objective',fun,'lb',[0 1e5 .5e5],'ub',[1 1e10 2e5],'x0',guess,'options',options);
%ms = MultiStart();
%[x,f] = run(ms,problem,100)

%------------------------------------------------------------------

function out=cost(A,B,C,D)
Anorm = A / trapz(div(:,1));
Bnorm = B / trapz(B(:,1));
Cnorm = C / trapz(nondiv(:,1));
Dnorm = D / trapz(D(:,1));

out1=1e5*immsre(Anorm,Bnorm);
out2=1e5*immsre(Cnorm,Dnorm);
out=out1+out2;

end

end