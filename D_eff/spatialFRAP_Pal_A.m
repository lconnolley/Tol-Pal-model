function [palA]=spatialFRAP_Pal_A(Dc,Db,alpha,~,gamma,kon,koff,N,sigma)

%P_d=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Fitting/Pal_dividing.mat');
A=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/TolPal/tolA_dividing.mat');

lngth=cellfun('size',A.cells,1);
L=median(lngth)*A.pixelsize;

%constants and parameters
x=-L/2:0.005*L:L/2;
t=A.t;

Df=Dc;              
Dp=0.000;
beta0_A=0;

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
bleach_d=interp1(-L/2:0.02*L:L/2,A.bleach,x);

%find steady state solution for ICs
[w1A,w2A,w3A,w4A]=steady_state_d(Dc,Db,alpha,beta0_A,gamma,kon,koff,N,sigma,L);

m=0;
sol = pdepe(m,@pdes,@ic,@bc,x,t);
cvA = sol(:,:,5);
fvA = sol(:,:,6);
pvA = sol(:,:,7);

palA=cvA+fvA+pvA;

disp('spatialFRAP_Pal_A completed')

% --------------------------------------------------------------

function [c,f,s] = pdes(x,~,w,DwDx)

c=[1; 1; 1; 1; 1; 1; 1];
f=[Dc; Db; Df; Dp; Dc; Df; Dp].*DwDx; 
s=[+alpha*w(2)*w(3) - beta0_A*beta(mu,x)*w(1) - gamma*w(1);   %Total complex
   -alpha*w(2)*w(3) + beta0_A*beta(mu,x)*w(1) + gamma*w(1);   %Total TolB
   -alpha*w(2)*w(3) + beta0_A*beta(mu,x)*w(1) + gamma*w(1) - kon*w(3)*(N-w(4)) + koff*w(4);   %Total free Pal
   +kon*w(3)*(N-w(4)) - koff*w(4);                        %Total bound Pal
   
   +alpha*w(2)*w(6) - beta0_A*beta(mu,x)*w(5) - gamma*w(5);   %Visible complex
   -alpha*w(2)*w(6) + beta0_A*beta(mu,x)*w(5) + gamma*w(5) - kon*w(6)*(N-w(4)) + koff*w(7);   %Visible free Pal
   +kon*w(6)*(N-w(4)) - koff*w(7)];                       %Visible bound Pal
end

% --------------------------------------------------------------

    function u = ic(x)
        
u=[w1A(end, find(xcopy==x,1));
    w2A(end, find(xcopy==x,1));
    w3A(end, find(xcopy==x,1));
    w4A(end, find(xcopy==x,1));
            
    w1A(end, find(xcopy==x,1)).*(bleach_d(find(xcopy==x,1)))';
    w3A(end, find(xcopy==x,1)).*(bleach_d(find(xcopy==x,1)))';
    w4A(end, find(xcopy==x,1)).*(bleach_d(find(xcopy==x,1)))'];

end

% --------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
    
pl = [0; 0; 0; 0; 0; 0; 0];     %negative=flux out, positive=flux in
ql = [1; 1; 1; 1; 1; 1; 1]; 
pr = pl;                        %negative=flux in, positive=flux out
qr = ql;

end

end