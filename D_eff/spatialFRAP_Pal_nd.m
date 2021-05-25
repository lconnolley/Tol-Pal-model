function [pal0]=spatialFRAP_Pal_nd(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma)

P=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/non-dividing.mat');

lngth=cellfun('length',P.cells);
L=mean(lngth)*P.pixelsize;

%constants and parameters
x=-L/2:0.005*L:L/2;
t=P.t;

Df=Dc;
Dp=0.000;

%define xcopy to get around limitations defining initial conditions
xcopy = x;

%extend bleach to match finer length scale
bleach_nd=interp1(-L/2:0.02*L:L/2,P.bleach,x);

%find steady state solution for ICs
[w10,w20,w30,w40]=steady_state_nd(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma,L);

m=0;
sol = pdepe(m,@pdes,@ic,@bc,x,t);
cv0 = sol(:,:,5);
fv0 = sol(:,:,6);
pv0 = sol(:,:,7);

pal0=cv0+fv0+pv0;

disp('spatialFRAP_Pal_nd completed')

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

u0=[w10(end, find(xcopy==x,1));
    w20(end, find(xcopy==x,1));
    w30(end, find(xcopy==x,1)); 
    w40(end, find(xcopy==x,1));

    
    w10(end, find(xcopy==x,1)).*(bleach_nd(find(xcopy==x,1)))';
    w30(end, find(xcopy==x,1)).*(bleach_nd(find(xcopy==x,1)))';
    w40(end, find(xcopy==x,1)).*(bleach_nd(find(xcopy==x,1)))'];

end

% --------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
    
pl = [0; 0; 0; 0; 0; 0; 0];     %negative=flux out, positive=flux in
ql = [1; 1; 1; 1; 1; 1; 1]; 
pr = pl;                        %negative=flux in, positive=flux out
qr = ql;

end

end