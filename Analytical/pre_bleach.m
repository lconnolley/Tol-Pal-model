function [c0,b0] = pre_bleach(d,a,b,t)

%--------------------------------------------------------------------------
%
%                            Complex, TolB
%
%--------------------------------------------------------------------------

%constants and parameters
L=1;
x=-L/2:0.01:L/2;


m=0;
sol0 = pdepe(m,@pdes0,@ic,@bc,x,t);      
c0 = sol0(:,:,1);                       
b0 = sol0(:,:,2);

C=c0(end,:);
B=b0(end,:);

save('steady_state.mat','C','B');

%{
%figures
figure(1)
subplot(2,1,1)
plot(x,c(end,:))
title('Distribution of TolB-Pal complex')
subplot(2,1,2)
plot(x,b(end,:))
title('Distribution of TolB in the IP')

%Amount of TolB in OP
div=trapz(x,c(end,:));
nondiv=trapz(x,c0(end,:));

%trapz(x,c(end,:)+b(end,:))
%}
% --------------------------------------------------------------

function [c,f,s] = pdes0(~,~,w,DwDx)
c=[1; 1];
f=[d; 1].*DwDx; 
s = [+a*w(2) - a*b*w(1);
     -a*w(2) + a*b*w(1)];
end

% --------------------------------------------------------------

function u0 = ic(~)
u0=[100; 5900];
end

% --------------------------------------------------------------

function [pl,ql,pr,qr] = bc(~,~,~,~,~)
pl = [0; 0];    %negative=flux out, positive=flux in
ql = [1; 1]; 
pr = pl;        %negative=flux in, positive=flux out
qr = ql;
end

end