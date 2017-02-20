%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              STEP 4b  - TOFTS                   %
%          Calculates Ktrans and Ve               %
%   using fminsearch to estimate residuals from   %
% convolution of the residue function and the AIF %
%                                                 %
%                                                 %
% Sara Maria Sprinkhuizen - March 2013         	  %
% based on code of                                %
% Kim Mouridsen  <kimm@nmr.mgh.harvard.edu>       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kin_par est_conc] = Step4b_Simplex_Fit_Kinetic_Tofts(obs_conc, AIF, init_params, FR_mins)


obs_conc;

t = 0:(length(obs_conc)-1);
t = t.*(FR_mins);

option=optimset('fminsearch');
option.MaxIter=1000000;
option.MaxFunEvals = 1000000;
% option.Display = 'iter'
kin_par = fminsearch(@kin_fit, init_params, option);
% kin_par = fminsearch(@kin_fit, init_params);

% k1=exp(kin_par(1));
% Ve=1/(1+exp(-kin_par(2)));
% k2 = k1/Ve;
% 
% res = exp(-k2.*t);
% est_conc = k1.*conv(AIF, res).*FR_mins;
% est_conc = est_conc(1:length(t));

% ktrans = exp(kin_par(1));
% ve = 1/(1+exp(-kin_par(2)));
% kep=ktrans/ve;
% t;
% est_conc=zeros(1,length(AIF));
% est_conc(1);
% for i=2:size(AIF,1)
%    log_e=-kep*t(2);
%    e=exp(log_e);
%    terma=AIF(i)*(e-log_e-1);
%    termb=AIF(i-1)*(e-(e*log_e)-1);
%    integral=(terma-termb)/power(log_e,2);
%    est_conc(i)=est_conc(i-1)*e + ktrans*t(2)*integral;
% end


  function cost = kin_fit(params)
    
    ktrans = exp(params(1));
    ve = 1/(1+exp(-params(2)));
    kep=ktrans/ve;
    t;
    est_conc=zeros(1,length(AIF));
    est_conc(1);
    for i=2:size(AIF,1)
       log_e=-kep*t(2);
       e=exp(log_e);
       terma=AIF(i)*(e-log_e-1);
       termb=AIF(i-1)*(e-(e*log_e)-1);
       integral=(terma-termb)/power(log_e,2);
       est_conc(i)=est_conc(i-1)*e + ktrans*t(2)*integral;
    end
   
      
%     k1=exp(params(1));
%     Ve=1/(1+exp(-params(2)));
%     k2 = k1/Ve;
% %     [k1, Ve, k2]
%     res = exp(-k2.*t);
%     est_conc = k1.*conv(AIF, res).*FR_mins;
%     est_conc = est_conc(1:length(t));

%     conv = polynomial convolution
%     est_conc
%     obs_conc

    if size(est_conc,2)~=1
        est_conc=est_conc';
    end
    d = (obs_conc-est_conc).^2;
    cost = sum(d);
    cost;
  end

end


