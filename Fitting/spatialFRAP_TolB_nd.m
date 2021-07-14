function tolb=spatialFRAP_TolB_nd(a,b,beta0)

B=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_nondiv_2s.mat');

lngth=cellfun('size',B.cells,1);
L=median(lngth)*B.pixelsize;

%find steady state solution for ICs
[w1,w2,w3,w4]=steady_state_nd(a,b,beta0,L);

disp('Starting TolB pdepe solver...')

%constants and parameters
lngth=cellfun('length',B.cells);
L=median(lngth)*B.pixelsize;
x=-L/2:0.005*L:L/2;
t=B.t;

Dc=(b*a)/(b-1);
Db=a/(b-1);
Df=Dc;
Dp=0.000;
alpha=5.4e-5;
gamma=0.006;
kon=1e-3;
koff=1;
N=1.9e5;

%define variable to get around limitations defining initial conditions
xcopy=x;

%extend bleach to match finer length scale
bleach=interp1(-L/2:0.02*L:L/2,B.bleach,x);

m=0;
%options=odeset('RelTol',100,'AbsTol',1);
sol=pdepe(m,@pdes,@ic,@bc,x,t);
cv = sol(:,:,5); %bleached complex
bv = sol(:,:,6); %bleached tolb

[m,~]=size(cv);
if m < length(t)
    tolb = ones(length(x),length(t)+1);
else
    factor=trapz(cv(1,:)+bv(1,:))/trapz(w1+w2);
    tolb=[(w1+w2)*factor; cv+bv]';
end


disp('TolB pdepe solver completed.')

%{
figure(3)
imagesc(t,x,tolb)
title('TolB')

figure(4)
subplot(2,1,1)
plot(x,tolb(:,1))
subplot(2,1,2)
plot(x,tolb(:,end))
%}

%--------------------------------------------------------------------------

function [c,f,s] = pdes(~,~,w,DwDx)
    
c=[1; 1; 1; 1; 1; 1];
f=[Dc; Db; Df; Dp; Dc; Db].*DwDx; 
s=[+alpha*w(2)*w(3) - (beta0/L)*w(1) - gamma*w(1);       %Total complex
   -alpha*w(2)*w(3) + (beta0/L)*w(1) + gamma*w(1);       %Total TolB
   -alpha*w(2)*w(3) + (beta0/L)*w(1) + gamma*w(1) - kon*w(3)*(N-w(4)) + koff*w(4);       %Total free Pal
   +kon*w(3)*(N-w(4)) - koff*w(4);                       %Total bound Pal
   
   +alpha*w(6)*w(3) - (beta0/L)*w(5) - gamma*w(5);       %Visible complex
   -alpha*w(6)*w(3) + (beta0/L)*w(5) + gamma*w(5)];      %Visible TolB
   
end

% -------------------------------------------------------------------------

function u0 = ic(x)

u0=[w1(end, find(xcopy==x,1,'first'));
    w2(end, find(xcopy==x,1,'first'));
    w3(end, find(xcopy==x,1,'first'));
    w4(end, find(xcopy==x,1,'first'));
    
    w1(end, find(xcopy==x,1,'first')).*(bleach(find(xcopy==x,1,'first')))';
    w2(end, find(xcopy==x,1,'first')).*(bleach(find(xcopy==x,1,'first')))'];

end

% -------------------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
    
pl = [0; 0; 0; 0; 0; 0];     %negative=flux out, positive=flux in
ql = [1; 1; 1; 1; 1; 1]; 
pr = pl;                     %negative=flux in, positive=flux out
qr = ql;

end

end