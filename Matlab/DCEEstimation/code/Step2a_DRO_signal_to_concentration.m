%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   STEP 2a                       %
%   Calculation of Contrast Agent Concentration   %
%     from the MRI Signal Intensity changes.	  %
%                                                 %
%                                                 %
% Sara Maria Sprinkhuizen - March 2013         	  %
% based on code of                                %
% Kim Mouridsen  <kimm@nmr.mgh.harvard.edu>       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gd_conc_4D, relSignal_4D] = Step2a_signal_to_concentration(input_signal_4D, T1_for_fit, firstBaseline, lastBaseline, r1_Cagent, TR, alpha_rad, smooth) 

signal_4D = input_signal_4D;

baselineVol_3D = mean(signal_4D(:,:,firstBaseline:lastBaseline),3);
relSignal_4D = zeros(size(signal_4D));
for i=1:size(signal_4D,3)
    relSignal_4D(:,:,i) = signal_4D(:,:,i)./baselineVol_3D;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The Eq used here for [Gd] -->
%
% [Gd] = -1/(r_1*TR) * ln (                 S(t)/S(t=0) * (TERM)   -  1                 )
%                         (-------------------------------------------------------------)
%                         (  exp(-TR*R1pre) * [ S(t)/S(t=0) * (TERM) * cos(alpha) -  1] )
% with:
%
% TERM = (      1 - exp(-TR*R1pre)        )
%        (--------------------------------)
%        ( 1 - cos(alpha)* exp(-TR*R1pre) )
%
%
%   This is derived from:  
%   S(t=postGd)        S_eq * sin(alpha) * (1 - exp{-TR(R1pre+r_1*[Gd])}) / (1 - cos(alpha)*exp{-TR(R1pre+r_1*[Gd])})
%   -----------   =     -------------------------------------------------------------------------------------------
%   S(t=preGd)              S_eq * sin(alpha) * (1 - exp{-TR*R1pre}) / (1 - cos(alpha)*exp{-TR*R1pre})
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R1pre = 1./T1_for_fit;   %1/msec
a = exp(-TR.*R1pre);
TERM = (1-a)./(1-a.*cos(alpha_rad));
y_4D = relSignal_4D.*(repmat(TERM,[size(signal_4D,1),size(signal_4D,2),size(signal_4D,3)])); 

% Use y_4D to calculate CA concentration:
gd_conc_4D = zeros(size(relSignal_4D));
for i=1:size(relSignal_4D,3);
    y_3D = squeeze(y_4D(:,:,i));
    gd_log_term = (y_3D-1)./(a.*(y_3D.*cos(alpha_rad)-1));
    gd_conc_4D(:,:,i) = -(1./(r1_Cagent*TR)) .* log(gd_log_term);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IS THE DRO DATA THE CONCENTRATION ALREADY ?!?!?!?
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gd_conc_4D = input_signal_4D;





