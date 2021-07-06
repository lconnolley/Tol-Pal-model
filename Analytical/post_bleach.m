function [c0,b0] = post_bleach(d,a,b,t)

%--------------------------------------------------------------------------
%
%                            Complex, TolB
%
%--------------------------------------------------------------------------

SS=load('steady_state.mat');

%constants and parameters
L=1;
x=-L/2:0.01:L/2;

%Photobleaching gaussian
sigma_f=0.1;
mu=0;
FRAP=@(mu,x) normpdf((x-mu)/sigma_f)/sigma_f/(normcdf((L-mu)/sigma_f)-normcdf(-mu/sigma_f));%truncated normal
p = max(FRAP(mu,x));
FRAP=@(mu,x) 1 - 0.5/p .* normpdf((x-mu)/sigma_f)/sigma_f/(normcdf((L-mu)/sigma_f)-normcdf(-mu/sigma_f));

%figure(1)
%clf
%plot(x,FRAP(mu,x))

xcopy = x;

m=0;
sol0 = pdepe(m,@pdes0,@ic,@bc,x,t);      
c0 = sol0(:,:,1);                       
b0 = sol0(:,:,2);

%figure(2)
%clf
%imagesc(c0')

trapz(c0)
% --------------------------------------------------------------

function [c,f,s] = pdes0(~,~,w,DwDx)
c=[1; 1];
f=[d; 1].*DwDx; 
s = [+a*w(2) - a*b*w(1);
     -a*w(2) + a*b*w(1)];
end

% --------------------------------------------------------------

function u0 = ic(x)
u0=[SS.C(end, find(xcopy==x,1,'first')).*FRAP(mu,x); 
    SS.B(end, find(xcopy==x,1,'first')).*FRAP(mu,x)];
end

% --------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
pl = [0; 0];    %negative=flux out, positive=flux in
ql = [1; 1]; 
pr = pl;        %negative=flux in, positive=flux out
qr = ql;
end

end