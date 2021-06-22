function pal = spatialFRAP_tolA_d(a,b,beta0)

A=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolA_dividing.mat');

lngth=cellfun('size',A.cells,1);
L=median(lngth)*A.pixelsize;

%find steady state solution for ICs
[w1,w2,w3,w4]=steady_state_tolA(a,b,0,L);

disp('Starting tolA pdepe solver...')

%constants and parameters
x=-L/2:0.005*L:L/2;
t=A.t; 


Dc=(b*a)/(b-1);
Db=a/(b-1);           
Df=Dc;
Dp=0.000;
alpha=5.4e4;
beta0=0;%tolA mutant, no sink
gamma=0.006;
kon=1e-3;
koff=1;       
N=1.7e5;

%define xcopy to get around limitations defining initial conditions
xcopy = x;

%extend bleach to match finer length scale
bleach=interp1(-L/2:0.02*L:L/2,A.bleach,x);

%shape of sink, beta
mu=0;
sigma=0.05;
beta=@(mu,x) normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma));%truncated normal
i = trapz(x,beta(mu,x));
beta=@(mu,x) 2/i*normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma)); %normalise to 2

q = trapz(x,beta(mu,x));
if q<1.99 || q>2.01
    error('Integral of beta function not equal to one.')
end
    
m=0;
%options=odeset('RelTol',1e-6,'AbsTol',1e-8);
sol = pdepe(m,@pdes,@ic,@bc,x,t);      
cv = sol(:,:,5);
fv = sol(:,:,6);
pv = sol(:,:,7);

[m,~]=size(cv);
if m < length(t)
    pal=ones(length(x),length(t)+1);
else
    factor=trapz(cv(1,:)+fv(1,:)+pv(1,:))/trapz(w1+w3+w4);
    pal=[(w1+w3+w4)*factor;cv+fv+pv]';
end

%{
figure(1)
clf
subplot(2,1,1)
plot(pal(:,1))
subplot(2,1,2)
plot(pal(:,end))

figure(2)
clf
imagesc(t,x,pal)
title('tolA')
%}


disp('tolA pdepe solver completed.')


% --------------------------------------------------------------

function [c,f,s] = pdes(x,~,w,DwDx)

c=[1; 1; 1; 1; 1; 1; 1];
f=[Dc; Db; Df; Dp; Dc; Df; Dp].*DwDx; 
s=[+alpha*w(2)*w(3) - beta0*beta(mu,x)*w(1) - gamma*w(1);   %Total complex
   -alpha*w(2)*w(3) + beta0*beta(mu,x)*w(1) + gamma*w(1);   %Total TolB
   -alpha*w(2)*w(3) + beta0*beta(mu,x)*w(1) + gamma*w(1) - kon*w(3)*(N-w(4)) + koff*w(4);   %Total free Pal
   +kon*w(3)*(N-w(4)) - koff*w(4);                        %Total bound Pal
   
   +alpha*w(2)*w(6) - beta0*beta(mu,x)*w(5) - gamma*w(5);   %Visible complex
   -alpha*w(2)*w(6) + beta0*beta(mu,x)*w(5) + gamma*w(5) - kon*w(6)*(N-w(4)) + koff*w(7);   %Visible free Pal
   +kon*w(6)*(N-w(4)) - koff*w(7)];                       %Visible bound Pal

end

% --------------------------------------------------------------

function u0 = ic(x)

u0=[w1(end, find(xcopy==x,1,'first'));
    w2(end, find(xcopy==x,1,'first'));
    w3(end, find(xcopy==x,1,'first')); 
    w4(end, find(xcopy==x,1,'first'));

    
    w1(end, find(xcopy==x,1,'first')).*(bleach(find(xcopy==x,1,'first')))';
    w3(end, find(xcopy==x,1,'first')).*(bleach(find(xcopy==x,1,'first')))';
    w4(end, find(xcopy==x,1,'first')).*(bleach(find(xcopy==x,1,'first')))'];

end

% --------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
    
pl = [0; 0; 0; 0; 0; 0; 0];     %negative=flux out, positive=flux in
ql = [1; 1; 1; 1; 1; 1; 1]; 
pr = pl;                        %negative=flux in, positive=flux out
qr = ql;

end

end