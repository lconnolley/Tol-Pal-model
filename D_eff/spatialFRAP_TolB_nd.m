function tolb=spatialFRAP_TolB_nd(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma)

B=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_nondiv_2s.mat');

lngth=cellfun('size',B.cells,1);
L=median(lngth)*B.pixelsize;

%constants and parameters
x=-L/2:0.005*L:L/2;
t=B.t;

Df=Dc;              %Victor's paper
Dp=0.00;

%define variable to get around limitations defining initial conditions
xcopy=x;

%extend bleach to match finer length scale
bleach=interp1(-L/2:0.02*L:L/2,B.bleach,x);

%find steady state solution for ICs
[w1,w2,w3,w4]=steady_state_nd(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma,L);

m=0;
sol=pdepe(m,@pdes,@ic,@bc,x,t);
cv = sol(:,:,5); %visible complex
bv = sol(:,:,6); %visible tolb

tolb=cv+bv;

%--------------------------------------------------------------------------

function [c,f,s] = pdes(~,~,w,DwDx)
    
c=[1; 1; 1; 1; 1; 1];
f=[Dc; Db; Df; Dp; Dc; Db].*DwDx; 
s=[+alpha*w(2)*w(3) - (beta0/L)*w(1) - gamma*w(1);       %Total complex
   -alpha*w(2)*w(3) + (beta0/L)*w(1) + gamma*w(1);       %Total TolB
   -alpha*w(2)*w(3) + (beta0/L)*w(1) + gamma*w(1) - kon*w(3)*(N-w(4)) + koff*w(4);       %Total free Pal
   +kon*w(3)*(N-w(4)) - koff*w(4);                            %Total bound Pal
   
   +alpha*w(6)*w(3) - (beta0/L)*w(5) - gamma*w(5);       %Visible complex
   -alpha*w(6)*w(3) + (beta0/L)*w(5) + gamma*w(5)];      %Visible TolB
   
end

%--------------------------------------------------------------------------

function u0 = ic(x)

u0=[w1(end, find(xcopy==x,1));
    w2(end, find(xcopy==x,1));
    w3(end, find(xcopy==x,1));
    w4(end, find(xcopy==x,1));
    
    w1(end, find(xcopy==x,1)).*(bleach(find(xcopy==x,1)))';
    w2(end, find(xcopy==x,1)).*(bleach(find(xcopy==x,1)))'];

end


%--------------------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
    
pl = [0; 0; 0; 0; 0; 0];     %negative=flux out, positive=flux in
ql = [1; 1; 1; 1; 1; 1]; 
pr = pl;                     %negative=flux in, positive=flux out
qr = ql;

end

end