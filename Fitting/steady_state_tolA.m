function [w1,w2,w3,w4] = steady_state_tolA(a,b,beta0,L)

%--------------------------------------------------------------------------
%
%     TolB, complex, Pal bound and unbound (kon//koff)   not M-M   
%
%--------------------------------------------------------------------------

%constants and parameters
x=-L/2:0.005*L:L/2;
t=0:60*2:60*20; 

disp('Finding steady state...')

Dc=(b*a)/(b-1);
Db=a/(b-1);   
Df=Dc;
Dp=0.000;
alpha=5.4e-5;
beta0=0;%tolA mutant, no sink
gamma=0.006;
kon=1e-4;
koff=1;
N=3.2e5;


%shape of sink, beta
mu=0;
sigma=0.08;
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
plot(x,w1+w2)

figure(2)
clf
plot(x,w1+w3+w4)

figure(3)
clf
imagesc(t,x,w1'+w3'+w4')
title('Pal')

figure(4)
clf
imagesc(t,x,w1'+w2')
title('TolB')
%}

figure(10)
x=-1/2:0.005:1/2;
clf
plot(x,w1,'DisplayName','TolB - outer')
hold on
plot(x,w2,'DisplayName','TolB - inner')
hold off
ylim([0 max(max(w1), max(w2))+500])
legend

disp('Steady state completed')

% --------------------------------------------------------------

function [c,f,s] = pdes(x,~,w,DwDx)

c=[1; 1; 1; 1];
f=[Dc; Db; Df; Dp].*DwDx; 
s=[+alpha*w(2)*w(3) - beta0*beta(mu,x)*w(1) - gamma*w(1);       %Complex
   -alpha*w(2)*w(3) + beta0*beta(mu,x)*w(1) + gamma*w(1);       %TolB
   -alpha*w(2)*w(3) + beta0*beta(mu,x)*w(1) + gamma*w(1) - kon*w(3)*(N - w(4)) + koff*w(4);   %Free Pal
   +kon*w(3)*(N - w(4)) - koff*w(4)];                            %Bound Pal

end

% --------------------------------------------------------------

function u0 = ic(~)

u0=[0; 9.96e3; 0; 99.6e3];

end

% --------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
pl = [0; 0; 0; 0];    %negative=flux out, positive=flux in
ql = [1; 1; 1; 1]; 
pr = pl;              %negative=flux in, positive=flux out
qr = ql;
end

end