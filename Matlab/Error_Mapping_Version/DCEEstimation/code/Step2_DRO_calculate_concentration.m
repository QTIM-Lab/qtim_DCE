%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   STEP 2                        %
%         TAILORED FOR THE DRO PROCESSING         %
%                                                 %
%                                                 %
%        Prepares for the calculation of          %
%         Contrast Agent Concentration            %
%      from the DRO Signal Intensity changes.	  %
%                                                 %
%                                                 %
% Sara Maria Sprinkhuizen - March 2013         	  %
% based on code of                                %d
% Kim Mouridsen  <kimm@nmr.mgh.harvard.edu>       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate concentration
function Step2_DRO_calculate_concentration(PreprocMapsCurrVisitDir,patID,MATfiles_dirname,TE1,TE2,TR,T1_fixed,alpha_rad,r1_Cagent,firstbaseline,lastbaseline,DCE_OPTS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT        
% - PreprocMapsCurrVisitDir : Folder containing prep files (dce_mc_*,T1,PD)
% - patID           : Example: 'CED_0101'
% - MATfiles_dirname: folder where CED_0X0X.mat file will be saved ind
% - TE1 and TE2     : echo times sequence millisec
% - TR              : repetition time sequence (is NOT the Frame Rate!) in millisec
% - T1_fixed        : should NOT be used here! Just in case of corrupt T1 map!
% - alpha_rad       : flipangle in radians, is fed to Step2a
% - r1_Cagent       : relaxivity of the Contrast agent. Hardcoded in Step0
% - firstbaseline   : first dynamic/frame volume to use for baseline signal calc
% - lastbaseline    : last dynamic/frame volume to use for baseline signal calc
% - DCE_OPTS        : User chosen processing options
%
% OUTPUT
% - Saves concentration volumes (conc_4D) and brainmask (amongst other things)
%   in a CED_0X0X.mat file in the MAT_FILES folder inside the RUN folder
% - no_prep_flag    : Indicates if NO prep files were found (=1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD T1
T1_for_fit = 1000; %assumed T1 in tissue       

cd (PreprocMapsCurrVisitDir)
fprintf('dir=%s\n', PreprocMapsCurrVisitDir);
SI_DRO_struct = load_untouch_nii('S0_500_6s_0s_sigma5.nii'); 

SI_DRO = squeeze(double(SI_DRO_struct.img));  % 
SI_DRO = permute(SI_DRO,[2,1,3]);


%%Calculate C(t) based on T2* (non)corrected signal
    disp ('   Calculating concentration....')
       [conc_4D, relSignal_4D] = Step2a_DRO_signal_to_concentration(SI_DRO,T1_for_fit,firstbaseline,lastbaseline,r1_Cagent,TR,alpha_rad,0);


figure,
for d = 1:61
imshow(conc_4D(:,:,d),[])
pause(.1)
end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % SARA: decide which of these things needed to be saved!
% ncurves = 15 * size(SI_E1_4D,3);  % Select on average 15 curves per slice
% maxvol = max(R2star_4D,[],4);
% nan_log = isnan(maxvol);
% maxvol(nan_log) = 0;
% maxvol(isinf(maxvol))=0;
% maxvol_OneD = maxvol(:);
% maxvol_sorted = sort(maxvol_OneD, 'descend');
% thres = maxvol_sorted(ncurves,1);
%vessel_mask = maxvol > thres;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear DCE_OPTS 
        
        disp (['Saving DRO.mat file'])   
        save([MATfiles_dirname patID '.mat'], '*')
        % Sara: the following volumes are saved in this CED0x0x.mat file:
        % - conc_4D
        % - SI_E1_4D
        % - SI_E1_R2star_corr_4D
        % - R2star_4D
        % - relSignal_4D
        % - wholebrainmask
        % - T1_for_fit



