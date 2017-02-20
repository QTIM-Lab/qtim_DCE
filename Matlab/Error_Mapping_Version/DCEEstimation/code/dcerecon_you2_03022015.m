
% to reconstruct ktrans and ve for a dce filename at voxel x,y,z

addpath /autofs/cluster/qtim/users/you2/DCEEstimation/code/
addpath /autofs/cluster/mi2b2/4/you2/Utilities/NIfTI_matlab_20140122/


% load dce1_mc.nii.gz
dce4D=load_untouch_nii('/autofs/cluster/qtim/users/you2/LongitudinalRegistration2/TMZ_03/dce0301/dce1_mc.nii.gz');
xsize=size(dce4D.img,1);
ysize=size(dce4D.img,2);
zsize=size(dce4D.img,3);
tsize=size(dce4D.img,4);
time=(1:tsize)*1645;  % ms

% load dependent parameters
DCE_OPTS = [  0    1    0    1   0   1   1   2    1    1   1    1    1    0 ];
T1_for_fit = 1000; %assumed T1 in tissue       
T1_fixed = 1000; %in milliseconds - NOTE THAT THIS MAY BE USED FOR THE CONCENTRATION CALCULATION INSTEAD OF T1 MAPS 
r1_Cagent = 0.0039; %Relaxometry of Contrast Agent used jkc

TR  = 5; % Repetition time in millisec! 
alpha=30;alpha_rad = (pi/180).*alpha;
FR = 0.5; %jkc
FR_mins = FR/60;    % in minutes 
Total_scan_time_mins = 11; %jkc
nr_of_frames =  Total_scan_time_mins/FR_mins+1; 
firstbaseline = 2; 
firstBaseline = firstbaseline;
lastbaseline = 35; 
lastBaseline = lastbaseline;
lpbs=lastbaseline+1;     %dynamic #



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for a point (61,28,11), take the 4D signal
%x=61;y=28;z=11;    % p1: a point with flat curve
%x=63; y=33; z=11;  % p2: another point with flat curve
%x=60; y=27; z=11;  % p3: a point with ascending curve
s=dce4D.img(x,y,z,:);   
input_signal_4D=zeros(1,1,tsize);
for n=1:tsize
    input_signal_4D(1,1,n)=s(1,1,1,n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert 4D signal to 4D concentration (reference to
% Step2a_DRO_signal_to_concentration.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal_4D = input_signal_4D;
baselineVol_3D = mean(signal_4D(:,:,firstBaseline:lastBaseline),3);
relSignal_4D = zeros(size(signal_4D));
for i=1:size(signal_4D,3)
    relSignal_4D(:,:,i) = signal_4D(:,:,i)./baselineVol_3D;
end
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
% in the above 
%figure; plot(squeeze(signal_4D)); title('observed signal');
%figure; plot(squeeze(gd_conc_4D)); title('observed concentration');
%figure; plot(squeeze(relSignal_4D)); title('normalized signal (normalized to the baseline)');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate aif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AIF]=generateAIF(tsize,FR,lpbs);
%figure; plot(squeeze(AIF)); title('AIF');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit ktrans and Ve (simplex)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obs_conc=squeeze(gd_conc_4D);
auc=trapz(time, obs_conc)/10^6;   % time has a unit of ms, obs_conc has a unit of mmol, so auc has a unit a mol*s
init_params=[-2, 0.1]; % initialization
[kin_par est_conc] = Step4b_Simplex_Fit_Kinetic_Tofts(obs_conc, AIF, init_params, FR_mins);
k1=exp(kin_par(1)); ktrans=k1;
Ve=1/(1+exp(-kin_par(2)));
k2 = k1/Ve;
fprintf('at (%d, %d, %d), Ve=%f, ktrans=%f\n', x, y, z, Ve, ktrans);
