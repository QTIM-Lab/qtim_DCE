%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   STEP 4                        %
%          Prepares for the modeling of           %
%    kinetic parameters from the DCE MRI data     %
%       Processes the results into STATS and      %
%       (overlay) maps and saves outcome          %
%                                                 %
%                                                 %
% Sara Maria Sprinkhuizen - March 2013         	  %
% based on code of                                %
% Kim Mouridsen  <kimm@nmr.mgh.harvard.edu>       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Step4_DRO_evaluateToftsModel(patID,MATfiles_dirname,PreprocDataDir,ROIpath_fname,DCE_PROC_MAPSfiles_dirname,DCE_PROC_STATSfiles_dirname,DCE_OPTS,FR_mins,lpbs,AIF,init_pars_fit,firstbaseline,lastbaseline)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT  
% - MATfiles_dirname  : To access the CED_0X0X.mat files containing
%                       CA concentration volumes
% - patID             : To load the correct CED_0101.mat filename
% - PreprocDataDir    : Where preprocessed data is (like ROI files)
% - ROIpath_fname     : specific path to ROI case the ROI masking is used
% - DCE_PROC_MAPSfiles_dirname   : to store MAPS files
% - DCE_PROC_STATSfiles_dirname  : to store STATS files
% - DCE_OPTS          : user chosen processing options
% - FR_mins           : Frame rate/dynamic scan time in minutes
% - lpbs              : Frame number of frame @ which CA bolus arrives
% - AIF               : Arterial Inpt Function
% - init_pars_fit     : Initial fitting parameters for Ktrans, Ve, (Vp) 
% - firstbaseline,lastbaseline: Used to calculate STDEV in baseline frames
%                       (before bolus) to do SMART mapping.
%
% OUTPUT
% - Processes and saves Ktrans, Ve, Kep, (Vp) data, maps, STATS
% - calcconcflag in case concentration calculation was not performed yet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(MATfiles_dirname)

    disp (['   Loading ' patID '.mat file'])   
    load(['DRO.mat']);
    % Sara: the following volumes are loaded from this CED0x0x.mat file:
    % - conc_4D (depending on the run, conc_4D is R2* corrected or not!)
    % - SI_E1_4D
    % - SI_E1_R2star_corr_4D
    % - R2star_4D
    % - relSignal_4D
    % - wholebrainmask
    % - T1_for_fit (not masked! do that here:) 
    
%%%%%%%% MASKING TYPE CHECK
MDL = 1;  %TOFTS
        
%%%%%%%%% Ktrans, Ve, Kep, estimated concentration curves are calculated for every voxel in certain masked region:
[Ktransmap, Ve, K2, Vp, est_conc_mat, obs_conc_mat, Curves_N_Pars_allmaskvoxels] = Step4a_DRO_Kinmap(conc_4D,FR_mins,lpbs,AIF,init_pars_fit,[],MDL,DCE_OPTS);
        
        
%%%%%%%%% Process output calculations:

% Concentration curves stats:
        STATS.est_conc.median = median(est_conc_mat);
        STATS.obs_conc.median = median(obs_conc_mat);
        STATS.est_conc.mean = mean(est_conc_mat);
        STATS.obs_conc.mean = mean(obs_conc_mat);
        % quality control jpgs:
            f=figure('Visible','off');plot(STATS.obs_conc.median,'b');hold on;plot(STATS.est_conc.median,'r'); 
            saveas(f,[DCE_PROC_STATSfiles_dirname patID '_mediancurve'],'jpg') 
            f=figure('Visible','off');plot(STATS.obs_conc.mean,'b');hold on;plot(STATS.est_conc.mean,'r')
            saveas(f,[DCE_PROC_STATSfiles_dirname patID '_meancurve'],'jpg') 
               
        disp(['Following data are for ' patID])
        disp ('=========================================================')
        disp(['Saving the Processed DCE output file for ' patID])
        disp ('=========================================================')
       

% Sara: save Ktrans, Ve, Kep, T1 only! Instead of all local variables, to save time and space
        fileName = [ DCE_PROC_MAPSfiles_dirname  patID '_KtransVeK2T1_MAPS.mat'];
        if MDL == 1 % Tofts
            save(fileName, 'Ktransmap', 'Ve', 'K2', 'T1_for_fit'); 

        elseif MDL == 2 % Extended Tofts
            save(fileName, 'Ktransmap', 'Ve', 'K2', 'Vp', 'T1_for_fit')
        end
        
% SARA: SAVE NIFTI
        Ktrans_nii = make_nii(Ktransmap);  
        fileName = [ DCE_PROC_MAPSfiles_dirname 'NII/'  patID '_Ktrans.nii']; save_nii(Ktrans_nii,fileName);
        Ve_nii = make_nii(Ve);
        fileName = [ DCE_PROC_MAPSfiles_dirname 'NII/' patID '_Ve.nii']; save_nii(Ve_nii,fileName);
        K2_nii = make_nii(K2);
        fileName = [ DCE_PROC_MAPSfiles_dirname 'NII/' patID '_K2.nii']; save_nii(K2_nii,fileName);
        T1_nii = make_nii(T1_for_fit);
        fileName = [ DCE_PROC_MAPSfiles_dirname 'NII/' patID '_T1.nii']; save_nii(T1_nii,fileName);


        
%% Sara added: save the sturct containing the observed and estimated concentration curves of ALL mask voxels
fileName      = [DCE_PROC_MAPSfiles_dirname(1:end-5) 'ALL_VOX/'  patID '_Curves_N_Pars_allmaskvoxels.mat'];
save(fileName, 'Curves_N_Pars_allmaskvoxels')

