function [Ct_ll]= estimateCtLogLinear(aif, FR, ktrans, ve)

numT=length(aif);
dt=ones(1,numT);
dt=FR*dt;
dt=dt/60;

kep=ktrans/ve;

Ct_ll=zeros(1,1200);
for i=2:size(dt,2)
   log_e=-kep*dt(i-1);
   e=exp(log_e);
   terma=aif(i)*(e-log_e-1);
   termb=aif(i-1)*(e-(e*log_e)-1);
   
   integral=(terma-termb)/power(log_e,2);
   Ct_ll(i)=Ct_ll(i-1)*e + ktrans*dt(i-1)*integral;
   
    
end

% figure, plot(Ct_ll)