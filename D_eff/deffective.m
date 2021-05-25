function [Deff]=deffective(t,data)

pixelsize=0.005;%size of each division
guess=1e-2;

[d]=fitkymo(t,data,guess); %fit
d=d*pixelsize^2;

Deff=d./data(:,1)/length(data(:,1));%scaled deff

end