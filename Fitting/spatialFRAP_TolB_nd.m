function tolb=spatialFRAP_TolB_nd(Dc,Db,beta0,N)

B=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_nondiv_2s_150s.mat');

%find steady state solution for ICs
[w1,w2,w3,w4]=steady_state_nd(Dc,Db,beta0,N);

%disp('Starting TolB pdepe solver...')

%constants and parameters
L=1.5;
x=-L/2:0.01*L:L/2;
t=B.t;

%Dc=0.02;
%Db=0.008;
Df=Dc;              %Victor's paper
Dp=0.000;
alpha=5.4e4;        %Papadakos paper
%beta0=9.48e8;
gamma=0.006;        %Papadakos paper
kon=1e-3;            %estimate
koff=1;            %estimate
%N=1.53e5;

%define variable to get around limitations defining initial conditions
xcopy=x;

%extend bleach to match finer length scale
bleach=interp1(-L/2:0.02*L:L/2,B.bleach,x);

m=0;
%options=odeset('RelTol',100,'AbsTol',1);
sol=pdepe(m,@pdes,@ic,@bc,x,t);
c = sol(:,:,1); %visible complex
b = sol(:,:,2); %visible tolb
f = sol(:,:,3);  %free pal
p = sol(:,:,4);  %bound pal

cv = sol(:,:,5); %bleached complex
bv = sol(:,:,6); %bleached tolb

[m,n]=size(c);
if m < length(t)
    tolb = ones(length(t),length(x));
else
    tolb=cv+bv;
end
tolb=tolb';


%disp('TolB pdepe solver completed.')

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
s=[+alpha*w(2)*w(3) - beta0*w(1) - gamma*w(1);       %Total complex
   -alpha*w(2)*w(3) + beta0*w(1) + gamma*w(1);       %Total TolB
   -alpha*w(2)*w(3) + beta0*w(1) + gamma*w(1) - kon*w(3)*(N-w(4)) + koff*w(4);       %Total free Pal
   +kon*w(3)*(N-w(4)) - koff*w(4);                            %Total bound Pal
   
   +alpha*w(6)*w(3) - beta0*w(5) - gamma*w(5);       %Visible complex
   -alpha*w(6)*w(3) + beta0*w(5) + gamma*w(5)];      %Visible TolB
   
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