function pal = spatialFRAP_Pal_nd(Dc,Db,beta0,N)

P=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/non-dividing.mat');

%find steady state solution for ICs
[w1,w2,w3,w4]=steady_state_nd(Dc,Db,beta0,N);

disp('Starting Pal pdepe solver...')

%constants and parameters
L=1.5;
x=-L/2:0.01*L:L/2;
t=P.t; 


%Dc=0.02;          %Dc>0.008
%Db=0.002;           
Df=Dc;              %Victor's paper
Dp=0.000;
alpha=5.4e4;        %Papadakos paper
%beta0=9.48e8;
gamma=0.006;        %Papadakos paper
kon=1e-3;            %Colin's estimate was 1e5-1e6
koff=1;            %Colin's estimate was 1-10
%N=1.53e5;

%define xcopy to get around limitations defining initial conditions
xcopy = x;

%extend bleach to match finer length scale
bleach=interp1(-L/2:0.02*L:L/2,P.bleach,x);
    
m=0;
%options=odeset('RelTol',1e-6,'AbsTol',1e-8);
sol = pdepe(m,@pdes,@ic,@bc,x,t);      
c = sol(:,:,1); %total complex
b = sol(:,:,2); %total tolb
f = sol(:,:,3); %total free pal
p = sol(:,:,4); %total bound pal

cv = sol(:,:,5); %visible complex
fv = sol(:,:,6); %visible free pal
pv = sol(:,:,7); %visible bound pal

[m,n]=size(c);
if m < length(t)
    pal=zeros(length(t),length(x));
else
    pal=cv+fv+pv;
end
pal=pal';

%{
%figure(1)
%subplot(2,1,1)
%plot(x,pal(:,1))
%subplot(2,1,2)
%plot(x,pal(:,end))

%figure(2)
%imagesc(t,x,pal)
%title('Pal')
%}

disp('Pal pdepe solver completed.')

% --------------------------------------------------------------

function [c,f,s] = pdes(~,~,w,DwDx)

c=[1; 1; 1; 1; 1; 1; 1];
f=[Dc; Db; Df; Dp; Dc; Df; Dp].*DwDx; 
s=[+alpha*w(2)*w(3) - beta0*w(1) - gamma*w(1);   %Total complex
   -alpha*w(2)*w(3) + beta0*w(1) + gamma*w(1);   %Total TolB
   -alpha*w(2)*w(3) + beta0*w(1) + gamma*w(1) - kon*w(3)*(N-w(4)) + koff*w(4);   %Total free Pal
   +kon*w(3)*(N-w(4)) - koff*w(4);                        %Total bound Pal
   
   +alpha*w(2)*w(6) - beta0*w(5) - gamma*w(5);   %Visible complex
   -alpha*w(2)*w(6) + beta0*w(5) + gamma*w(5) - kon*w(6)*(N-w(4)) + koff*w(7);   %Visible free Pal
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