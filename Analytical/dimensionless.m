clear all

x=-1/2:0.01:1/2;

[Bin_s, Bout_s] = conc_profile_s(50,100,1/5);
[Bin_ns, Bout_ns] = conc_profile_ns(50,100,1/5);

figure(1)
clf
plot(x,Bin_s,'color',[0, 0.447, 0.741],'DisplayName','Bin - sink')
hold on
plot(x,Bout_s,'color',[0.85, 0.325, 0.098],'DisplayName','Bout - sink')
yline(Bin_ns,'--','color',[0, 0.447, 0.741],'DisplayName','Bin - no sink')
yline(Bout_ns,'--','color',[0.85, 0.325, 0.098],'DisplayName','Bout - no sink')
hold off
legend

d=trapz(x,Bout_s)
nd=Bout_ns

%%
%Varying d

n=[];
s=[];
for d=0.01:0.01:10

    [NS,S] = tot_conc(50,100,1/d);
    n=[n;NS];
    s=[s;S];
    
end

figure(2)
clf
plot(0.01:0.01:10,n,'DisplayName','No sink','color',[0.8500, 0.3250, 0.0980])
hold on
plot(0.01:0.01:10,s,'DisplayName','Sink','color',[0 0.4470 0.7410])
hold off
legend
title('Varying the ratio of the diffusion rates, d')

%%
%Varying b

n=[];
s=[];
for b=0:0.1:150

    [NS,S] = tot_conc(50,b,1/5);
    n=[n;NS];
    s=[s;S];
    
end

figure(3)
clf
plot(0:0.1:150,n,'DisplayName','No sink','color',[0.8500, 0.3250, 0.0980])
hold on
plot(0:0.1:150,s,'DisplayName','Sink','color',[0 0.4470 0.7410])
hold off
legend
title('Varying the rate of exchange, b')

%%
%Varying a

n=[];
s=[];
for a=0:1:300

    [NS,S] = tot_conc(a,100,1/5);
    n=[n;NS];
    s=[s;S];
    
end

figure(4)
clf
plot(0:1:300,n,'DisplayName','No sink','color',[0.8500, 0.3250, 0.0980])
hold on
plot(0:1:300,s,'DisplayName','Sink','color',[0 0.4470 0.7410])
hold off
legend
title('Varying the rate of exchange,a')

%%
%varying b and d

dd=[(0.001:0.2:2.001)'; (2:0.5:5)']';
aa=[(0:0.5:5)'; (5:5:40)']';

n=[];
s=[];
for d=dd
    for b=aa
        
        [N,S]=tot_conc(50,b,1/d);
        n=[n;N];
        s=[s;S];
        
    end
end

bsz=length(aa);%block size
l=length(n);%total length of data

cell_n=mat2cell(n,diff([0:bsz:l-1,l]));
cell_s=mat2cell(s,diff([0:bsz:l-1,l]));

cell_n=[cell_n{:}];
cell_n=cell_n';
cell_s=[cell_s{:}];
cell_s=cell_s';

figure(5)
clf
X=aa;
Y=dd;
Z1=cell_s-cell_n;
Z2=cell_s;
Z3=cell_n;
surf(X,Y,Z2)
hold on
surf(X,Y,Z3)
hold off
view(142,19)
ax=gca;
ax.FontSize=20;
xticks([0 10 20 30 40])
yticks([0 1 2 3 4 5])
%ylabel('Ratio of diffusion rates, d')
xlabel('Transport rate, b')
zlabel('Total concentration of Bout')

%%
%varying a and d
%{
dd=0:0.2:1.5;
aa=0:25:600;

n=[];
s=[];
for d=dd
    for a=aa
        
        [N,S]=tot_conc(a,10,d);
        n=[n;N];
        s=[s;S];
        
    end
end

bsz=length(aa);%block size
l=length(n);%total length of data

cell_n=mat2cell(n,diff([0:bsz:l-1,l]));
cell_s=mat2cell(s,diff([0:bsz:l-1,l]));

cell_n=[cell_n{:}];
cell_n=cell_n';
cell_s=[cell_s{:}];
cell_s=cell_s';

figure(6)
clf
X=aa;
Y=dd;
Z1=cell_s-cell_n;
Z2=cell_s;
Z3=cell_n;
surf(X,Y,Z2)
hold on
surf(X,Y,Z3)
hold off
view(149,28)
ax=gca;
ax.FontSize=20;
ylabel('Ratio of diffusion rates, d')
xlabel('Transport rate, a')
zlabel('Total concentration of Aout')

%}
%%

function [Bin_s, Bout_s] = conc_profile_s(a,b,d)
%concentration profile in the inner and outer periplasm for a point 
%sink

x=-1/2:0.01:1/2;
kappa=sqrt(a);
kappabar=(kappa/2 * (1+cosh(kappa)) / (sinh(kappa))) - 1;

Bbar= 1 / (b + 1 + (a*b*kappabar/(d*kappa^2)));

G=(kappa/2)*(cosh(kappa*x) + cosh(kappa*(abs(x)-1))) / (sinh(kappa));

D=-d/(a*b) - 1/(2*kappa) * (1+cosh(kappa))/sinh(kappa);

Bin_s = b*Bbar*G;
Bout_s = - (a*b*Bbar)/d * (G/kappa^2 + D);

end

function [Bin_ns, Bout_ns] = conc_profile_ns(~,b,~)

Bin_ns = 1 / (1 + 1/b);
Bout_ns = 1 / (1 + b);

end

function [ns, s] = tot_conc(a,b,d)
%total concentration in the outer periplasm for no sink and a point sink

kappa=sqrt(a);
kappabar=(kappa/2)*((1+cosh(kappa))/(sinh(kappa))) - 1;

%calculations are for Bin, BT - Bin = Bout

ns= 1 - 1/(1+1/b);
s = 1 - 1 / (1 + 1/b + kappabar/d);

end