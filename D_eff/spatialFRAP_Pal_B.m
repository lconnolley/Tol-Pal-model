function [palB]=spatialFRAP_Pal_B(Dc,Db,~,beta0,gamma,kon,koff,N,sigma)

P=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Fitting/Pal_dividing.mat');
B=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolB_dividing.mat');

lngth=cellfun('length',B.cells);
L=mean(lngth)*B.pixelsize;

%constants and parameters
x=-L/2:0.005*L:L/2;
t=B.t;

Df=Dc;
Dp=0.000;
alpha_B=0;

%shape of sink, beta
mu=0;
beta=@(mu,x) normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma));%truncated normal
i = trapz(x,beta(mu,x));
beta=@(mu,x) 2/i*normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma)); %normalise to 2

q = trapz(x,beta(mu,x));
if q<1.99 || q>2.01
    error('Integral of beta function not equal to one.')
end

%define xcopy to get around limitations defining initial conditions
xcopy = x;

%extend bleach to match finer length scale
bleach_d=interp1(-L/2:0.02*L:L/2,P.bleach,x);

%find steady state solution for ICs
[w1B,w2B,w3B,w4B]=steady_state_d(Dc,Db,alpha_B,beta0,gamma,kon,koff,N,sigma,L);

m=0;
sol = pdepe(m,@pdes,@ic,@bc,x,t);      
cvB = sol(:,:,5);
fvB = sol(:,:,6);
pvB = sol(:,:,7);

palB=cvB+fvB+pvB;

disp('spatialFRAP_Pal_B completed')

% --------------------------------------------------------------

function [c,f,s] = pdes(x,~,w,DwDx)

c=[1; 1; 1; 1; 1; 1; 1];
f=[Dc; Db; Df; Dp; Dc; Df; Dp].*DwDx; 
s=[+alpha_B*w(2)*w(3) - beta0*beta(mu,x)*w(1) - gamma*w(1);   %Total complex
   -alpha_B*w(2)*w(3) + beta0*beta(mu,x)*w(1) + gamma*w(1);   %Total TolB
   -alpha_B*w(2)*w(3) + beta0*beta(mu,x)*w(1) + gamma*w(1) - kon*w(3)*(N-w(4)) + koff*w(4);   %Total free Pal
   +kon*w(3)*(N-w(4)) - koff*w(4);                        %Total bound Pal
   
   +alpha_B*w(2)*w(6) - beta0*beta(mu,x)*w(5) - gamma*w(5);   %Visible complex
   -alpha_B*w(2)*w(6) + beta0*beta(mu,x)*w(5) + gamma*w(5) - kon*w(6)*(N-w(4)) + koff*w(7);   %Visible free Pal
   +kon*w(6)*(N-w(4)) - koff*w(7)];                       %Visible bound Pal

end

% --------------------------------------------------------------
function u = ic(x)
        
u=[w1B(end, find(xcopy==x,1));
    w2B(end, find(xcopy==x,1));
    w3B(end, find(xcopy==x,1));
    w4B(end, find(xcopy==x,1));
            
    w1B(end, find(xcopy==x,1)).*(bleach_d(find(xcopy==x,1)))';
    w3B(end, find(xcopy==x,1)).*(bleach_d(find(xcopy==x,1)))';
    w4B(end, find(xcopy==x,1)).*(bleach_d(find(xcopy==x,1)))'];

end

% --------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
    
pl = [0; 0; 0; 0; 0; 0; 0];     %negative=flux out, positive=flux in
ql = [1; 1; 1; 1; 1; 1; 1]; 
pr = pl;                        %negative=flux in, positive=flux out
qr = ql;

end

end