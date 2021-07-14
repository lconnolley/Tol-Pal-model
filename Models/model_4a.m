 function model_4a()

%--------------------------------------------------------------------------
%
%     TolB, complex, Pal bound and unbound (kon//koff)   not M-M   
%
%--------------------------------------------------------------------------

%d=[0.0027, 6.656e-4, 8.28e8, 5.3e4, 1.28e5];%Dc, Db, beta0, N, kon
%d=[0.0027, 6.656e-4, 1.0e9, 5.3e4, 1.28e5];%Dc, Db, beta0, N, kon
%d=[0.0027, 0.0014, 7.95e8, 2e5, 1.0747e-3];%Dc, Db, beta0, N, kon
%d=[0.0027, 0.0014, 7.95e8, 3e6, 1e5];
d=[0.0192, 0.0057, 2.17e-2, 2e5, 1e-3];

%constants and parameters
L=3;    %length of an average e. coli cell
x_step=0.005;
x=-L/2:x_step*3:L/2;
x0=-L/4:x_step*L/2:L/4;
x1=-1/2:x_step:1/2;
t=0:10*60:60*60; 
Dc=d(1);            %Dc>0.0125
Db=d(2);           %Db<0.004
Df=Dc;              %Victor's paper
Dp=0.00;
alpha=5.4e-5;        %Papadakos paper
beta0=d(3);
gamma=0.006;        %Papadakos paper
kon=d(5);             %estimate was 1e5-1e6 M^-1 s^-1
koff=10;            %estimate was 1-10 s^-1
N=d(4);

%mutants
alpha_B=0;
beta0_P=0;

%shape of sink, beta
mu=0;
sigma=0.04*L;
beta=@(mu,x) normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma));%truncated normal
i = trapz(x,beta(mu,x));
beta=@(mu,x) 2/i*normpdf((x-mu)/sigma)/sigma/(normcdf((L-mu)/sigma)-normcdf(-mu/sigma)); %normalise to 2

%double check
q = trapz(x,beta(mu,x));
if q<1.99 || q>2.01
    error('Integral of beta function not equal to two.')
end
    
m=0;
sol = pdepe(m,@pdes,@ic,@bc,x,t);      
c = sol(:,:,1); 
b = sol(:,:,2);
f = sol(:,:,3);
p = sol(:,:,4);
disp('Dividing done.')

sol0 = pdepe(m,@pdes0,@ic,@bc,x0,t);      
c0 = sol0(:,:,1);                       
b0 = sol0(:,:,2);
f0 = sol0(:,:,3);
p0 = sol0(:,:,4);
disp('Non-dividing done.')

%{
sol_B = pdepe(m,@pdes_B,@ic,@bc,x0,t);      
c_B = sol_B(:,:,1);                       
b_B = sol_B(:,:,2);
f_B = sol_B(:,:,3);
p_B = sol_B(:,:,4);

sol_P = pdepe(m,@pdes_P,@ic,@bc,x0,t);      
c_P = sol_P(:,:,1);                       
b_P = sol_P(:,:,2);
f_P = sol_P(:,:,3);
p_P = sol_P(:,:,4);
%}
%figures
figure(4)
clf
subplot(4,1,1)
plot(x1,c0(end,:))
hold on
plot(x1,c(end,:))
hold off
title('Distribution of TolB-Pal complex')
subplot(4,1,2)
plot(x1,b0(end,:))
hold on
plot(x1,b(end,:))
hold off
title('Distribution of TolB in the IP')
subplot(4,1,3)
plot(x1,f0(end,:))
hold on
plot(x1,f(end,:))
hold off
title('Distribution of unbound Pal')
subplot(4,1,4)
plot(x1,p0(end,:))
hold on
plot(x1,p(end,:))
hold off
title('Distribution of bound Pal')

figure(5)
clf
subplot(2,1,1)
plot(x1,c0(end,:)+f0(end,:)+p0(end,:))
hold on
plot(x1,c(end,:)+f(end,:)+p(end,:))
hold off
title('Distribution of total Pal')
subplot(2,1,2)
plot(x1,c0(end,:)+b0(end,:))
hold on
plot(x1,c(end,:)+b(end,:))
hold off
title('Distribution of total TolB')

figure(6)
clf
plot(x1,(c0(end,:)*Dc + f0(end,:)*Df) ./ (c0(end,:)+f0(end,:)+p0(end,:)),'DisplayName','Non-dividing')
hold on
plot(x1,(c(end,:)*Dc + f(end,:)*Df) ./ (c(end,:)+f(end,:)+p(end,:)),'DisplayName','Dividing')
%plot(x1,(c_B(end,:)*Dc + f_B(end,:)*Df) ./ (c_B(end,:)+f_B(end,:)+p_B(end,:)),'DisplayName','TolB mutant')
%plot(x1,(c_P(end,:)*Dc + f_P(end,:)*Df) ./ (c_P(end,:)+f_P(end,:)+p_P(end,:)),'DisplayName','TolA mutant')
hold off
%title('Pal - Effective diffusion coefficient')
%ylim([0 1.2e-3])
xlabel('Relative position')
ylabel('Effective diffusion coefficient')
legend

figure(7)
clf
plot(x1,(c0(end,:)*Dc + b0(end,:)*Db) ./ (c0(end,:) + b0(end,:)))
hold on
plot(x1,(c(end,:)*Dc + b(end,:)*Db) ./ (c(end,:) + b(end,:)))
hold off
title('TolB - Effective diffusion coefficient')

%Amount of TolB in OP/IP in non-div
i=trapz(x,b(end,:));
o=trapz(x,c(end,:));

tolb_ip=i/(i+o)*100
tolb_op=o/(i+o)*100


% --------------------------------------------------------------

function [c,f,s] = pdes(x,t,w,DwDx)
c=[1; 1; 1; 1];
f=[Dc; Db; Df; Dp].*DwDx; 
s=[+(alpha*w(2)*w(3)) - (beta0*beta(mu,x))*w(1) - gamma*w(1);       %Complex
   -(alpha*w(2)*w(3)) + (beta0*beta(mu,x))*w(1) + gamma*w(1);       %TolB
   -(alpha*w(2)*w(3)) + (beta0*beta(mu,x))*w(1) + gamma*w(1) - kon*w(3)*(N - w(4)) + koff*w(4);   %Free Pal
   +kon*w(3)*(N - w(4)) - koff*w(4)]                               %Bound Pal
end

function [c,f,s] = pdes0(x0,t,w,DwDx)
c=[1; 1; 1; 1];
f=[Dc; Db; Df; Dp].*DwDx; 
s=[+(alpha*w(2)*w(3)) - (beta0/L)*w(1) - gamma*w(1);
   -(alpha*w(2)*w(3)) + (beta0/L)*w(1) + gamma*w(1);
   -(alpha*w(2)*w(3)) + (beta0/L)*w(1) + gamma*w(1) - kon*w(3)*(N - w(4)) + koff*w(4);
   +kon*w(3)*(N - w(4)) - koff*w(4)];
end

function [c,f,s] = pdes_B(x0,t,w,DwDx)
c=[1; 1; 1; 1];
f=[Dc; Db; Df; Dp].*DwDx; 
s=[+(alpha_B*w(2)*w(3)) - (beta0/L)*w(1) - gamma*w(1);
   -(alpha_B*w(2)*w(3)) + (beta0/L)*w(1) + gamma*w(1);
   -(alpha_B*w(2)*w(3)) + (beta0/L)*w(1) + gamma*w(1) - kon*w(3)*(N - w(4)) + koff*w(4);
   +kon*w(3)*(N - w(4)) - koff*w(4)];
end

function [c,f,s] = pdes_P(x0,t,w,DwDx)
c=[1; 1; 1; 1];
f=[Dc; Db; Df; Dp].*DwDx; 
s=[+(alpha*w(2)*w(3)) - (beta0_P/L)*w(1) - gamma*w(1);
   -(alpha*w(2)*w(3)) + (beta0_P/L)*w(1) + gamma*w(1);
   -(alpha*w(2)*w(3)) + (beta0_P/L)*w(1) + gamma*w(1) - kon*w(3)*(N - w(4)) + koff*w(4);
   +kon*w(3)*(N - w(4)) - koff*w(4)];
end

% --------------------------------------------------------------

function u0 = ic(~)
    
u0=[0; 9.96e3; 0; 99.6e3];

end

% --------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
pl = [0; 0; 0; 0];    %negative=flux out, positive=flux in
ql = [1; 1; 1; 1]; 
pr = pl;              %negative=flux in, positive=flux out
qr = ql;
end

end