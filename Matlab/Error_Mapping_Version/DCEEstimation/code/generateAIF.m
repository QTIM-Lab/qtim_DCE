function [AIF]=generateAIF(n,FR,lpbs)


%% Generate AIF
% t is the number of time points
% model choices are 1. Parker (default), 2. Fluckinger 3. FDA

%% times

timeOfBolus=lpbs;
AIF=zeros(n,1);
timePoints=n-timeOfBolus; % was 60-timeOfBolus. Yangming changed on 05/20/2014
%t=0.1*[0:timePoints-1];
t=FR*[0:timePoints-1]/60;

%% Parker
% defining parameters
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

%%
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

% 
% 
% %% Fluckinger
% 
% %% defining parameters
% 
% a1=6;
% delta1=1.0;
% alpha=2.92;
% tau1=0.0442;
% a2=1.1208;
% delta2=delta1+0.2227;
% tau2=0.1430;
% a3=0.3024;
% delta3=delta1+0.6083;
% a4=0.7164;
% T=7.8940;

%% 
