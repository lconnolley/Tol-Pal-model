function [d,fval]=fitkymo6(tolb_d,pal_d,tolb_nd,pal_nd,tolA,tolB,guess)

fun=@(d) cost(tolb_d(:,2:end), spatialFRAP_TolB_d(d(1),d(2),d(3),d(4)), pal_d(:,2:end), spatialFRAP_Pal_d(d(1),d(2),d(3),d(4)), tolb_nd(:,2:end), spatialFRAP_TolB_nd(d(1),d(2),d(3),d(4)), pal_nd(:,2:end), spatialFRAP_Pal_nd(d(1),d(2),d(3),d(4)), tolA(:,2:end), spatialFRAP_tolA_d(d(1),d(2),d(3),d(4)), tolB(:,2:end), spatialFRAP_Pal_d(d(1),d(2),0,d(4)));
%Pattern search
options=optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,'InitialMeshSize',0.1,'Display','iter');
[d,fval,exitflag,output] = patternsearch(fun,guess,[],[],[],[],[0 0 0 0],[],[],options);

%Surrogate optimisation
%options = optimoptions('surrogateopt','UseParallel',true,'MaxFunctionEvaluations',1000);
%[d,fval,exitflag,output] = surrogateopt(fun,[0 0 0 0],[10 10 1e20 1e20],options);

%Multistart
%rng default %for reproducibility
%options = optimoptions(@fmincon,'UseParallel',true,'Display','iter');
%problem = createOptimProblem('fmincon','objective',fun,'lb',[0 0 1e5 .5e5],'ub',[1 1 1e10 2e5],'x0',guess,'options',options);
%ms = MultiStart();
%[x,f] = run(ms,problem,100)

%------------------------------------------------------------------

function out=cost(A,B,C,D,E,F,G,H,I,J,K,L)
Anorm = A / trapz(tolb_d(:,1));
Bnorm = B / trapz(B(:,1));
Cnorm = C / trapz(pal_d(:,1));
Dnorm = D / trapz(D(:,1));
Enorm = E / trapz(tolb_nd(:,1));
Fnorm = F / trapz(F(:,1));
Gnorm = G / trapz(pal_nd(:,1));
Hnorm = H / trapz(H(:,1));
Inorm = I / trapz(tolA(:,1));
Jnorm = J / trapz(J(:,1));
Knorm = K / trapz(tolB(:,1));
Lnorm = L / trapz(L(:,1));

out1=1e5*immsre(Anorm,Bnorm);
out2=1e5*immsre(Cnorm,Dnorm);
out3=1e5*immsre(Enorm,Fnorm);
out4=1e5*immsre(Gnorm,Hnorm);
out5=1e5*immsre(Inorm,Jnorm);
out6=1e5*immsre(Knorm,Lnorm);
out=out1+out2+out3+out4+out5+out6;

end

end