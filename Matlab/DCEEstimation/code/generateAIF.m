function [AIF]=generateAIF(tsize,time_interval_seconds,bolus_arrival_time)

% This function uses the Parker model to generate a population-based AIF
% for a certain duration and bolus arrival time.

timeOfBolus=bolus_arrival_time;
AIF=zeros(tsize,1);
timePoints=tsize-timeOfBolus;
t=time_interval_seconds*[0:timePoints-1]/60;

% Parker model parameters
% TODO: Add source.
a1=0.809;
a2=0.330;
T1=0.17406;
T2=0.365;
sigma1=0.0563;
sigma2=0.132;
alpha=1.050;
beta=0.1685;
s=38.078;
tau=0.483;

% Parker model calculations.
term0=alpha*exp(-beta*t)./(1+exp(-s*(t-tau)));
term1=[];
term2=[];
A1=a1/(sigma1*((2*pi)^0.5));
B1=exp(-(t-T1).^2./(2.*sigma1^2));
term1=A1.*B1;
A2=a2/(sigma2*((2*pi)^0.5));
B2=exp(-(t-T2).^2./(2.*sigma2^2));
term2=A2.*B2;
aifPost=term0+term1+term2;
sp=timeOfBolus+1;
AIF(sp:end)=aifPost;
