function [Deff,d]=deffective(t,data)
%finds the best fit to the SpatialFRAP data

guess=1e-2;

[d]=fitkymo(t,data,guess); %fit

Deff=d./data(:,1)/length(data(:,1));%scaled deff

end