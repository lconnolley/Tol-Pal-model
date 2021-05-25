function tolb=spatialFRAP_TolB(~)

B_d=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_dividing_5s.mat');
B_nd=load('/home/connolleyl/Documents/ownCloud/Tol-Pal/MATLAB/Import/TolB_nondividing_5s.mat');

%constants and parameters
L=3;
x_step=0.01;
x=-L/2:x_step*L:L/2;
x0=-L/4:x_step*L/2:L/4;
x1=-1/2:x_step:1/2;
t=0:10:200; 

Dc=0.0102;          %Dc>0.0125
Db=0.0056;          %Db<0.004
Df=Dc;              %Victor's paper
Dp=0.00;
alpha=5.4e4;        %Papadakos paper
beta0=1.28e9;
gamma=0.006;        %Papadakos paper
kon=1e3;            %estimate was 1e5-1e6
koff=10;            %estimate was 1-10
N=1.62e5;

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

%define variable to get around limitations defining initial conditions
xcopy=x;

%extend bleach to match finer length scale
bleach_d=interp1(-L/2:0.02*L:L/2,B_d.bleach,x);
bleach_nd=interp1(-L/2:0.02*L:L/2,B_nd.bleach,x);

%find steady state solution for ICs
[w1,w2,w3,w4]=steady_state_d(Dc,Db,alpha,beta0,gamma,kon,koff,N,sigma);
[w10,w20,w30,w40]=steady_state_nd(Dc,Db,alpha,beta0,gamma,kon,koff,N);

m=0;
sol=pdepe(m,@pdes,@ic,@bc,x,t);
%c = sol(:,:,1); %total complex
%b = sol(:,:,2); %total tolb
%f = sol(:,:,3); %total free pal
%p = sol(:,:,4); %total bound pal
cv = sol(:,:,5); %visible complex
bv = sol(:,:,6); %visible tolb

sol0=pdepe(m,@pdes0,@ic0,@bc,x,t);
cv0 = sol0(:,:,5); %visible complex
bv0 = sol0(:,:,6); %visible tolb

tolb=cv+bv;
tolb0=cv0+bv0;

figure(20)
clf
subplot(3,1,1)
plot(x1,w10+w20,'DisplayName','Non-dividing')
hold on
plot(x1,w1+w2,'DisplayName','Dividing')
hold off
title('TolB pre bleach')
subplot(3,1,2)
plot(x1,tolb0(1,:),'DisplayName','Non-dividing')
hold on
plot(x1,tolb(1,:),'DisplayName','Dividing')
hold off
title('TolB after bleaching')
legend
subplot(3,1,3)
plot(x1,tolb0(end,:),'DisplayName','Non-dividing')
hold on
plot(x1,tolb(end,:),'DisplayName','Dividing')
hold off
title('TolB after recovery')
legend


data=[w1+w2; tolb]';
data0=[w10+w20; tolb0]';

Deff=deffective(t,data);
Deff0=deffective(t,data0);

figure(21)
clf
plot(x1,Deff0,'DisplayName','Non-dividing')
hold on
plot(x1,Deff,'DisplayName','Dividing')
hold off
xlabel('Relative position')
ylabel('Effective diffusion coefficient')
legend
ylim([0 0.045])
%--------------------------------------------------------------------------

function [c,f,s] = pdes(x,~,w,DwDx)
    
c=[1; 1; 1; 1; 1; 1];
f=[Dc; Db; Df; Dp; Dc; Db].*DwDx; 
s=[+alpha*w(2)*w(3) - beta0*beta(mu,x)*w(1) - gamma*w(1);       %Total complex
   -alpha*w(2)*w(3) + beta0*beta(mu,x)*w(1) + gamma*w(1);       %Total TolB
   -alpha*w(2)*w(3) + beta0*beta(mu,x)*w(1) + gamma*w(1) - kon*w(3)*(N-w(4)) + koff*w(4);       %Total free Pal
   +kon*w(3)*(N-w(4)) - koff*w(4);                            %Total bound Pal
   
   +alpha*w(6)*w(3) - beta0*beta(mu,x)*w(5) - gamma*w(5);       %Visible complex
   -alpha*w(6)*w(3) + beta0*beta(mu,x)*w(5) + gamma*w(5)];      %Visible TolB
   
end

function [c,f,s] = pdes0(~,~,w,DwDx)
    
c=[1; 1; 1; 1; 1; 1];
f=[Dc; Db; Df; Dp; Dc; Db].*DwDx; 
s=[+alpha*w(2)*w(3) - beta0*w(1) - gamma*w(1);       %Total complex
   -alpha*w(2)*w(3) + beta0*w(1) + gamma*w(1);       %Total TolB
   -alpha*w(2)*w(3) + beta0*w(1) + gamma*w(1) - kon*w(3)*(N-w(4)) + koff*w(4);       %Total free Pal
   +kon*w(3)*(N-w(4)) - koff*w(4);                            %Total bound Pal
   
   +alpha*w(6)*w(3) - beta0*w(5) - gamma*w(5);       %Visible complex
   -alpha*w(6)*w(3) + beta0*w(5) + gamma*w(5)];      %Visible TolB
   
end

%--------------------------------------------------------------------------

function u0 = ic(x)

u0=[w1(end, find(xcopy==x,1));
    w2(end, find(xcopy==x,1));
    w3(end, find(xcopy==x,1));
    w4(end, find(xcopy==x,1));
    
    w1(end, find(xcopy==x,1)).*(bleach_nd(find(xcopy==x,1)))';
    w2(end, find(xcopy==x,1)).*(bleach_nd(find(xcopy==x,1)))'];

end

function u0 = ic0(x)

u0=[w10(end, find(xcopy==x,1));
    w20(end, find(xcopy==x,1));
    w30(end, find(xcopy==x,1));
    w40(end, find(xcopy==x,1));
    
    w10(end, find(xcopy==x,1)).*(bleach_d(find(xcopy==x,1)))';
    w20(end, find(xcopy==x,1)).*(bleach_d(find(xcopy==x,1)))'];

end

%--------------------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
    
pl = [0; 0; 0; 0; 0; 0];     %negative=flux out, positive=flux in
ql = [1; 1; 1; 1; 1; 1]; 
pr = pl;                     %negative=flux in, positive=flux out
qr = ql;

end

end
