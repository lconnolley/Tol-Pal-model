function [A1, A2, A10, A20] = steadystate_model2(d,a,b)

%constants and parameters
x=-1/2:0.01:1/2;
t=0:5:60;

%d
%a
%b

%shape of sink, beta
mu=0;
sigma=0.05;
beta=@(mu,x) normpdf((x-mu)/sigma)/sigma/(normcdf((1-mu)/sigma)-normcdf(-mu/sigma));%truncated normal
i = trapz(x,beta(mu,x));
beta=@(mu,x) 1/i*normpdf((x-mu)/sigma)/sigma/(normcdf((1-mu)/sigma)-normcdf(-mu/sigma)); %normalise to 1

q = trapz(x,beta(mu,x));
if q<0.99 || q>1.01
    error('Integral of beta function not equal to one.')
end

m=0;
sol = pdepe(m,@pdes,@ic,@bc,x,t);      
A1 = sol(:,:,1); 
A2 = sol(:,:,2);

sol0 = pdepe(m,@pdes0,@ic,@bc,x,t);      
A10 = sol0(:,:,1);                       
A20 = sol0(:,:,2);

%{
figure(1)
clf
imagesc(t,x,A1')

figure(2)
clf
imagesc(t,x,A2')
%}

%steady state solutions
A1 = A1(end,:);
A2 = A2(end,:);

A10 = A10(end,:);
A20 = A20(end,:);


% --------------------------------------------------------------

function [c,f,s] = pdes(x,~,w,DwDx)

c = [1; 1];
f = [d; 1].*DwDx; 
s = [+a*w(2) - b*beta(mu,x)*w(1);
     -a*w(2) + b*beta(mu,x)*w(1)];
end

function [c,f,s] = pdes0(~,~,w,DwDx)
c=[1; 1];
f=[d; 1].*DwDx; 
s = [+a*w(2) - b*w(1);
     -a*w(2) + b*w(1)];
end

% --------------------------------------------------------------

function u0 = ic(~)
    
u0=[0; 1];

end

% --------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
pl = [0; 0];    %negative=flux out, positive=flux in
ql = [1; 1]; 
pr = pl;        %negative=flux in, positive=flux out
qr = ql;
end

end