function sol=spatialFRAP(Dc)

SS=load('steady_state.mat');
P=load('Pal_dividing.mat');

%constants and parameters
x=0:0.02:1;
t=0:2*60:10*60;

disp('Starting pdepe solver...')

%Dc=0.0125;
Db=0.002;
Df=Dc;              %Victor's paper
Dp=0.00;
alpha=5.4e4;        %Papadakos paper
beta0=4e8;
gamma=0.006;        %Papadakos paper
kon=1e3;            %estimate
koff=10;            %estimate
N=2e5;

%options=odeset('RelTol',1e-6,'AbsTol',1e-12);
sol = pdepe(0,@pdes,@ic,@bc,x,t);
sol1=sol(:,:,1);
sol2=sol(:,:,2);
sol3=sol(:,:,3);

sol=sol1+sol2+sol3;
sol=sol';

disp('pdepe solver completed.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------
function [c,f,s] = pdes(xx,t,w,dudx)
c = [1; 1; 1];
f= [Dc; Db; 0].*dudx;
s = [+alpha*w(2)*w(3) - beta0*w(1) - gamma*w(1);
    -alpha*w(2)*w(3) + beta0*w(1) + gamma*w(1);
    -alpha*w(2)*w(3) + beta0*w(1) + gamma*w(1)];
end
% --------------------------------------------------------------
function out = ic(xx)
out = [0; 6000; 60000];
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bc(xl,ul,xr,ur,t)
pl = [0; 0; 0]; %negative=flux out, positive= flux in
ql = [1; 1; 1]; 
pr = pl; % negative=flux in, positive=flux out
qr = ql;
end


end