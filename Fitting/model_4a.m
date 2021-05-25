function model_4a()

%--------------------------------------------------------------------------
%
%     TolB, complex, Pal bound and unbound (kon//koff)   not M-M   
%
%--------------------------------------------------------------------------

%constants and parameters
L=1;
x=-L/2:0.02:L/2;
t=0:10*60:60*60; 

Dc=0.009; %Dc=0.0125;
Db=0.005; %Db=0.01;
Df=Dc;              %Victor's paper
Dp=0.00;
alpha=5.4e4;        %Papadakos paper
beta0=2e7;
gamma=0.006;        %Papadakos paper
kon=1e4;            %Colin's estimate was 1e5-1e6
koff=10;            %Colin's estimate was 1-10
N=1.7e5;


%shape of sink, beta
mu=0;
sigma=0.03;
beta=@(mu,x) 1/2*normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma));%truncated normal

q = trapz(x,beta(mu,x));
if q<0.99 || q>1.01
    error('Integral of beta function not equal to one.')
end
    
m=0;
sol = pdepe(m,@pdes,@ic,@bc,x,t);      
c = sol(:,:,1); 
b = sol(:,:,2);
f = sol(:,:,3);
p = sol(:,:,4);

sol0 = pdepe(m,@pdes0,@ic,@bc,x,t);      
c0 = sol0(:,:,1);                       
b0 = sol0(:,:,2);
f0 = sol0(:,:,3);
p0 = sol0(:,:,4);

%Steady state solutions
w1 = c(end,:);
w2 = b(end,:);
w3 = f(end,:);
w4 = p(end,:);
save('steady_state','w*')


%figures
figure(7)
clf
subplot(4,1,1)
plot(x,c(end,:))
title('Distribution of TolB-Pal complex')
subplot(4,1,2)
plot(x,b(end,:))
title('Distribution of TolB in the IP')
subplot(4,1,3)
plot(x,f(end,:))
title('Distribution of unbound Pal')
subplot(4,1,4)
plot(x,p(end,:))
hold on
plot(x,p0(end,:))
hold off
title('Distribution of bound Pal')

figure(8)
clf
subplot(2,1,1)
plot(x,c(end,:)+f(end,:)+p(end,:))
hold on
plot(x,c0(end,:)+f0(end,:)+p0(end,:))
hold off
title('Distribution of total Pal')
subplot(2,1,2)
plot(x,c(end,:)+b(end,:))
title('Distribution of total TolB')

figure(9)
clf
plot(x,(c(end,:)*Dc + f(end,:)*Df) ./ (c(end,:)+f(end,:)+p(end,:)))
hold on
plot(x,(c0(end,:)*Dc + f0(end,:)*Df) ./ (c0(end,:)+f0(end,:)+p0(end,:)))
hold off
title('Effective diffusion coefficient')
%ylim([0 1.2e-3])


%{
figure(4)
pcolor(x,t,p+c+f)
shading interp

%Amount of TolB in OP
tolb_div=trapz(x,c(end,:))
tolb_nondiv=trapz(x,c0(end,:))

%Deff mean
deff_div=mean((c(end,:)*Dc) ./ (c(end,:)+p(end,:)+A(end,:)))
deff_nondiv=mean((c0(end,:)*Dc) ./ (c0(end,:)+p0(end,:)+A(end,:)))

TolB = trapz(x,c(end,:)+b(end,:))
%}

% --------------------------------------------------------------

function [c,f,s] = pdes(x,t,w,DwDx)

c=[1; 1; 1; 1];
f=[Dc; Db; Df; Dp].*DwDx; 
s=[+(alpha*w(2)*w(3)) - beta0*beta(mu,x)*w(1) - gamma*w(1);       %Complex
   -(alpha*w(2)*w(3)) + beta0*beta(mu,x)*w(1) + gamma*w(1);       %TolB
   -(alpha*w(2)*w(3)) + beta0*beta(mu,x)*w(1) + gamma*w(1) - kon*w(3)*(1 - w(4)/N) + koff*w(4);   %Free Pal
   +kon*w(3)*(1 - w(4)/N) - koff*w(4)];                             %Bound Pal

 end

function [c,f,s] = pdes0(x,t,w,DwDx)

c=[1; 1; 1; 1];
f=[Dc; Db; Df; Dp].*DwDx; 
s=[+(alpha*w(2)*w(3)) - beta0*w(1) - gamma*w(1);       %Complex
   -(alpha*w(2)*w(3)) + beta0*w(1) + gamma*w(1);       %TolB
   -(alpha*w(2)*w(3)) + beta0*w(1) + gamma*w(1) - kon*w(3)*(1 - w(4)/N) + koff*w(4);   %Free Pal
   +kon*w(3)*(1 - w(4)/N) - koff*w(4)];                             %Bound Pal

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