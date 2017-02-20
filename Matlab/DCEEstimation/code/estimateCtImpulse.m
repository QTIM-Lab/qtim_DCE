function [Ct_imp]= estimateCtImpulse(aif, FR, ktrans, ve)

numT=length(aif);
dt=ones(1,numT);
dt=FR*dt;
dt=dt/60;
deltaT=FR/60;

kep=ktrans/ve;

e=exp(-kep.*dt);

Ct_imp(1)=0;
Ct_imp=zeros(1,1200);
for i=1:size(dt,2)-1
%     Ct_imp(i+1)=Ct_imp(i).*e(i)+kep*(dt(i).*aif(i+1));
      Ct_imp(i+1)=Ct_imp(i).*e(i)+ktrans*(deltaT.*aif(i+1));

    
end

% figure, plot(Ct_imp)
