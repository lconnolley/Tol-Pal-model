function [Bin,Bout] = non_div(a,b,d)

%variables
L=1;
x=-L/2:0.01:L/2;
t=0:60:60*60;


m=0;
sol=pdepe(m,@pdes,@ic,@bc,x,t);
Bin=sol(:,:,1);
Bout=sol(:,:,2);

time=size(t)
vec=size(Bout)

Bin=Bin(end,:);
Bout=Bout(end,:);

figure(5)
clf
plot(x,Bin,'DisplayName','Bin-pdepe')
hold on
plot(x,Bout,'DisplayName','Bout-pdepe')
yline(1/(1+1/b),'--','color',[0, 0.447, 0.741],'DisplayName','Bin')
yline(1/(1+b),'--','color',[0.85, 0.325, 0.098],'DisplayName','Bout')
hold off
legend

%{
figure(1)
clf
plot(x,Bin,'DisplayName','Bin')
hold on
plot(x,Bout,'DisplayName','Bout')
hold off
legend
%}
% -------------------------------------------------------------------------

    function [c,f,s] = pdes(~,~,w,DwDx)
        c=[1; 1];
        f=[1; d].*DwDx;
        s=[-a*w(1) + a*b*w(2); %Bin
           +a*w(1) - a*b*w(2)]; %Bout
    end

    function u0 = ic(~)
        u0 = [0; 1];
    end

    function [pl,ql,pr,qr] = bc(~,~,~,~,~)
        pl=[0; 0];
        ql=[1; 1];
        pr=pl;
        qr=ql;
    end

end