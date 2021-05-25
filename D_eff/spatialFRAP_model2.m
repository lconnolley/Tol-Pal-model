function nothing = spatialFRAP_model2()

%constants and parameters
x=-1/2:0.01:1/2;
q1=0.001;%0.005;
e1=0.005;%0.05;
t=0:q1:e1;

d=10;
a=80000;
b=40000;%320000;
a0=0;
b0=0;

%find steady state solution for ICs
[a1,a2,a10,a20]=steadystate_model2(d,a,b);
%[a1M,a2M,a10M,a20M]=steadystate_model2(d,a,b);

%define xcopy to get around limitations defining initial conditions
xcopy=x;

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

%shape of bleach
tau=0.4;
bleach=@(mu,x) normpdf((x-mu)/tau)/tau/(normcdf((1-mu)/sigma)-normcdf(-mu/tau));%truncated normal
i = max(bleach(mu,x));
bleach=@(mu,x) 1/i*normpdf((x-mu)/tau)/tau/(normcdf((1-mu)/tau)-normcdf(-mu/tau)); %normalise to 1
factor=1-trapz(x,bleach(mu,x));

%q = max(bleach(mu,x));
%if q<0.49 || q>.51
%   error('Integral of bleach function not equal to one.')
%end

m=0;
sol = pdepe(m,@pdes,@ic,@bc,x,t);      
A1 = sol(:,:,3);
A2 = sol(:,:,4);

sol0 = pdepe(m,@pdes0,@ic0,@bc,x,t);
A10 = sol0(:,:,3);
A20 = sol0(:,:,4);

solM = pdepe(m,@pdesM,@ic0,@bc,x,t);
A1M = solM(:,:,3);
A2M = solM(:,:,4);


figure(4)
clf
subplot(2,1,1)
plot(x,a1+a2)
subplot(2,1,2)
plot(x,A1(1,:)+A2(1,:))

figure(1)
clf

subplot(3,1,1)
imagesc(-q1:q1:e1,x,[(a1+a2)'.*factor',(A1+A2)'])
xlabel('Relative time')
ylabel('Relative position')
title('A1+A2 - point sink')

subplot(3,1,2)
imagesc(-q1:q1:e1,x,[(a10+a20)'.*factor',(A10+A20)'])
xlabel('Relative time')
ylabel('Relative position')
title('A1+A2 - no sink')

subplot(3,1,3)
imagesc(-q1:q1:e1,x,[(a10+a20)'.*factor',(A1M+A2M)'])  
xlabel('Relative time')
ylabel('Relative position')
title('a,b = 0')


data = [(a1+a2)'.*factor',(A1+A2)'];
data0 = [(a10+a20)'.*factor',(A10+A20)'];
dataM = [(a10+a20)'.*factor',(A1M+A2M)'];

Deff = deffective(t,data);
Deff0 = deffective(t,data0);
DeffM = deffective(t,dataM);

x1=-1/2:0.01:1/2;

figure(2)
clf
plot(x1,Deff,'DisplayName','Point sink')
hold on
plot(x1,Deff0,'DisplayName','No sink')
plot(x1,DeffM,'DisplayName','a,b=0')
%plot(x1,DeffM0,'DisplayName','No sink; a,b=0')
hold off
xlabel('Relative position')
legend;


Deff0(1)
DeffM(1)

%{
figure(5)
clf
plot(x,(A1(end,:).*d + A2(end,:)) ./ (A1(end,:)+A2(end,:)),'DisplayName','Point sink')
hold on
plot(x,(A10(end,:).*d + A20(end,:)) ./ (A10(end,:)+A20(end,:)),'DisplayName','No sink')
plot(x,(A1M(end,:).*d + A2M(end,:)) ./ (A1M(end,:)+A2M(end,:)),'DisplayName','a,b=0')
hold off
xlabel('Relative position')
ylabel('Weighted average diffusion constant')
legend;
%}
figure(3)
clf
imagesc(-q1:q1:e1,x,(A1+A2)') 

figure(6)
clf
plot(x,A1(4,:)+A2(4,:))

% --------------------------------------------------------------

function [c,f,s] = pdes(x,~,w,DwDx)

c=[1; 1; 1; 1];
f=[d; 1; d; 1].*DwDx;
s = [+a*w(2) - b*beta(mu,x)*w(1);      %total A1
     -a*w(2) + b*beta(mu,x)*w(1);      %total A2
     
     +a*w(4) - b*beta(mu,x)*w(3);      %visible A1
     -a*w(4) + b*beta(mu,x)*w(3)];     %visible A2
end

function [c,f,s] = pdes0(~,~,w,DwDx)

c=[1; 1; 1; 1];
f=[d; 1; d; 1].*DwDx;
s = [+a*w(2) - b*w(1);      %total A1
     -a*w(2) + b*w(1);      %total A2
     
     +a*w(4) - b*w(3);      %visible A1
     -a*w(4) + b*w(3)];     %visible A2
end

function [c,f,s] = pdesM(~,~,w,DwDx)

c=[1; 1; 1; 1];
f=[d; 1; d; 1].*DwDx;
s = [+a0*w(2) - b0*w(1);      %total A1
     -a0*w(2) + b0*w(1);      %total A2
     
     +a0*w(4) - b0*w(3);      %visible A1
     -a0*w(4) + b0*w(3)];     %visible A2
end

% --------------------------------------------------------------

function u0 = ic(x)

u0=[a1(end,find(xcopy==x,1,'first'));
    a2(end,find(xcopy==x,1,'first'));
    
    a1(end,find(xcopy==x,1,'first')).*(1-bleach(mu,x));
    a2(end,find(xcopy==x,1,'first')).*(1-bleach(mu,x))];

end

function u0 = ic0(x)

u0=[a10(end,find(xcopy==x,1,'first'));
    a20(end,find(xcopy==x,1,'first'));
    
    a10(end,find(xcopy==x,1,'first')).*(1-bleach(mu,x));
    a20(end,find(xcopy==x,1,'first')).*(1-bleach(mu,x))];
    
end

%{
function u0 = ic0M(x)

u0=[a10M(end,find(xcopy==x,1,'first'));
    a20M(end,find(xcopy==x,1,'first'));
    
    a10M(end,find(xcopy==x,1,'first')).*(1-bleach(mu,x));
    a20M(end,find(xcopy==x,1,'first')).*(1-bleach(mu,x))];
    
end
%}
% --------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
pl = [0; 0; 0; 0];    %negative=flux out, positive=flux in
ql = [1; 1; 1; 1]; 
pr = pl;              %negative=flux in, positive=flux out
qr = ql;
end

end