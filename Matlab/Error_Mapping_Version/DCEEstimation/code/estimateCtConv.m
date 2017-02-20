function [Ct_conv]= estimateCtConv(aif, FR, ktrans, ve,hem)

numT=length(aif);
FR_min=FR/60;
% t = ones(1, numT); 
t= 0:(numT-1);
% t = t.*FR/60;
t = t*FR_min;
kep=ktrans/ve;

res = exp(-kep*t);
Ct_conv = (1/(1-hem))*ktrans*conv(aif', res)*FR_min;% (1/(1-hem))*
Ct_conv = Ct_conv(1:length(t));


% figure, plot(Ct_conv)

