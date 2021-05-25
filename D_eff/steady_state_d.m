function [w1,w2,w3,w4] = steady_state_d(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma,L)

%--------------------------------------------------------------------------
%
%     TolB, complex, Pal bound and unbound (kon//koff)   not M-M   
%
%--------------------------------------------------------------------------

%constants and parameters
x=-L/2:0.005*L:L/2;
t=0:60*2:60*20; 
 
Df=Dc;
Dp=0.000;


%shape of sink, beta
mu=0;
%sigma=0.05;
beta=@(mu,x) normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma));%truncated normal
i = trapz(x,beta(mu,x));
beta=@(mu,x) 2/i*normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma)); %normalise to 2

q = trapz(x,beta(mu,x));
if q<1.99 || q>2.01
    error('Integral of beta function not equal to two.')
end

m=0;
sol = pdepe(m,@pdes,@ic,@bc,x,t);
c = sol(:,:,1);
b = sol(:,:,2);
f = sol(:,:,3);
p = sol(:,:,4);

%Steady state solutions
w1 = c(end,:);
w2 = b(end,:);
w3 = f(end,:);
w4 = p(end,:);

%{
figure(1)
clf
plot(w1+w3+w4)
%}
% --------------------------------------------------------------

function [c,f,s] = pdes(x,t,w,DwDx)

c=[1; 1; 1; 1];
f=[Dc; Db; Df; Dp].*DwDx; 
s=[+(alpha*w(2)*w(3)) - beta0*beta(mu,x)*w(1) - gamma*w(1);       %Complex
   -(alpha*w(2)*w(3)) + beta0*beta(mu,x)*w(1) + gamma*w(1);       %TolB
   -(alpha*w(2)*w(3)) + beta0*beta(mu,x)*w(1) + gamma*w(1) - kon*w(3)*(N - w(4)) + koff*w(4);   %Free Pal
   +kon*w(3)*(N - w(4)) - koff*w(4)];                             %Bound Pal

end

% --------------------------------------------------------------

function u0 = ic(x)

u0=[0; 9.96e3; 0; 99.6e3];

end

% --------------------------------------------------------------

function [pl,ql,pr,qr] = bc(xl,ul,xr,ur,t)
pl = [0; 0; 0; 0];    %negative=flux out, positive=flux in
ql = [1; 1; 1; 1]; 
pr = pl;              %negative=flux in, positive=flux out
qr = ql;
end

end