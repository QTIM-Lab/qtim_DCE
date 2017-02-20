function [AIF]=Step3_DRO_generateAIF(n,FR_mins,lpbs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT         
% - n             :  Number of Frames for final AIF = nr of dynamics
% - FR            :  Frame rate/Dynamic scan time: NEEDS TO BE MINUTES!
% - lpbs          :  Time of Bolus in Number of Frames
%
% OUTPUT
% - AIF           :  Arterial Input Function based on model of Parker 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% cd /autofs/cluster/qtim/users/sara/DCE_DRO/RUN_selS_selV_TE1_R2starcor_popAIF_Tfts_inpHC_subT1_FIG_ROI/MAT_FILES
cd /Users/kalpathy/Documents/MATLAB/dceDRO/RUN_selS_selV_TE1_R2starcor_popAIF_Tfts_inpHC_subT1_FIG_ROI/MAT_FILES
%% Generate AIF

load('DRO.mat')

AIFregion = SI_DRO(70:80,:,:);
AIF = squeeze(mean(mean(AIFregion,1),2));
figure,plot(AIF);title('AIF')