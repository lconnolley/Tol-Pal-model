function sol=spatialFRAP(t,y0,D)
x=1:1:length(y0);
xint=1;       %x(2)-x(1);
dDdx=gradient(D,1);


%options=odeset('RelTol',1e-6,'AbsTol',1e-12);
sol = pdepe(0,@pdes,@ic,@bc,x,t);
sol=sol';
% figure(5)
% imagesc(t,x,sol)
% 
% figure(6)
% plot(x,D)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------
function [c,f,s] = pdes(xx,t,u,dudx)
c = 1;
if length(D)==1
    f = D*dudx;
else
    index=round(xx/xint);
    f = dDdx(index)*u+D(index)*dudx;
end
s = 0;
end
% --------------------------------------------------------------
function out = ic(xx)
out = y0(xx)';
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = bc(~,~,~,~,~)
pl = [0]; %negative=flux out, positive= flux in
ql = [1]; 
pr = [0]; % negative=flux in, positive=flux out
qr = [1];
end


end