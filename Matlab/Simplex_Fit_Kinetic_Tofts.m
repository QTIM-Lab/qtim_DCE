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
% updated by Andrew Beers <abeers@mgh.harvard.edu>%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kinetic_params, estimated_concentration] = Simplex_Fit_Kinetic_Tofts(observed_concentration, gd_AIF, initial_params, FR_mins)

% Create a time series for the given data. The Tofts model assumes time is
% counted in minutes.
t = 0:(length(observed_concentration)-1);
t = t.*(FR_mins);

% fminsearch uses the Nelder-Mead a.k.a. Simplex a.k.a. Amoeba fitting
% algorithm to fit curves derived from the two-parameter Tofts model to 
% provided 'observed' curves. 
option=optimset('fminsearch');
option.MaxIter=10000;
option.MaxFunEvals = 10000;
kinetic_params = fminsearch(@kinetic_fit, initial_params, option);

% TODO: Verbose mode for debugging.
% option.Display = 'iter'

function cost = kinetic_fit(params)
    
    ktrans = exp(params(1));
    ve = 1/(1+exp(-params(2)));
    kep=ktrans/ve;

    estimated_concentration=zeros(1,length(gd_AIF));

    % Integration method 1. Calculated for each time-step.
    for i=2:size(gd_AIF,1)
       log_e=-kep*t(2);
       e=exp(log_e);
       terma=gd_AIF(i)*(e-log_e-1);
       termb=gd_AIF(i-1)*(e-(e*log_e)-1);
       integral=(terma-termb)/power(log_e,2);
       estimated_concentration(i)=estimated_concentration(i-1)*e + ktrans*t(2)*integral;
    end

    % Integration method 2. Calculated via convolution. Experimentally, has
    % provided less accurate results.
    % res = exp(-kep.*t);
    % estimated_concentration = ktrans.*conv(AIF, res).*FR_mins;
    % estimated_concentration = estimated_concentration(1:length(t));

    % Sanity check for bad input data.
    if size(estimated_concentration,2)~=1
        estimated_concentration=estimated_concentration';
    end
    
    % Cost function - squared error.
    d = (observed_concentration-estimated_concentration).^2;
    cost = sum(d);
end

end


