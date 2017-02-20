%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                 %
%          DCE MRI PROCESSING SOFTWARE            %
%       Tailored for the DRO-data                 %
%                                                 %
%         DIGITAL REFERENCE OBJECT                %
%                                                 %
% https://dblab.duhs.duke.edu/modules/            %  
%            QIBAcontent/index.php?id=1           %
%                                                 %
%                                                 %
% Four main Steps:                                %
%   1. Get Subject specific scan info/pars        %      
%   2. Calculate Concentration volumes            %
%   3. Calculate AIF                              %
%   4. Perform DCE kinetic modeling               %
%                                                 %
%                                                 %
%  Sara Maria Sprinkhuizen                        %
%  March 2013                                     %
%                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc

addpath /autofs/cluster/mi2b2/4/you2/Utilities/NIfTI_20140122

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT

%%%  INPUT: PATHS
% Set paths to patient data, preprocessed data and output folder
% - Select Directory With all Subject Data
%SubjectDir = '/autofs/cluster/qtim/users/sara/DCE_DRO';
SubjectDir = '/autofs/cluster/qtim/users/you2/DCEEstimation/test/'
%SubjectDir = '/Users/kalpathy/Documents/MATLAB/dceDRO/';
%PreprocDataDir  = '/autofs/cluster/qtim/users/sara/DCE_DRO';
PreprocDataDir  = '/autofs/cluster/qtim/users/you2/DCEEstimation/test/'
% PreprocDataDir  = '/Users/kalpathy/Documents/MATLAB/dceDRO/';

% - Where to save the output?
% outputVisitDir='/autofs/cluster/qtim/users/sara/DCE_DRO/';
% outputVisitDir='/Users/kalpathy/Documents/MATLAB/dceDRO/';
outputVisitDir='../test/';

%%%  INPUT: MANUAL POSTPROCESSING CHOICES 
% Manual input for the following choices/manual inputs:
% DCE_OPTS(1) = AllSubjects; 1 is yes, 0 is no
% DCE_OPTS(2) = SelectSubjects; 1 is yes, 0 is no
% DCE_OPTS(3) = AllVisits;
% DCE_OPTS(4) = SelectVisits;
% DCE_OPTS(5) = Baseline;
% DCE_OPTS(6) => 1 = TE1; 2 = TE2 --> NB (R2* corr SI TE1) = (R2* corr SI TE2)!
% DCE_OPTS(7) => 1 = R2* corrected; 0 = not R2* corrected
% DCE_OPTS(8) => 1 = Whole brain; 2 = T1AxialROI; 3 = SmartMap
% DCE_OPTS(9) => 1 = PopulationAIF; 2 = Subject specific AIF
% DCE_OPTS(10) => 1 = Tofts; 2 = Extended Tofts;
% DCE_OPTS(11) => 1 = Use hardcoded input_pars ; 2 = Use SVD method for input pars
% DCE_OPTS(12) => 1 = Use Subj specific T1 map; 2 = use fixed hardcoded T1
%                 3 = Use DCE-baseline-derived T1map 
% DCE_OPTS(13) = show graphics during process
% DCE_OPTS(14) => 0 = do not smooth SI before concentration calc, 1 = smooth before
 
%            allS selS allV selV bV TE  R2* Mask AIF  MDL inp   T1  FIG  smth];
%            (1)  (2)  (3)  (4) (5) (6) (7) (8)  (9) (10) (11) (12) (13) (14) 
DCE_OPTS = [  0    1    0    1   0   1   1   2    1    1   1    1    1    0 ];

%%% INPUT: HARDCODED INPUT
T1_fixed = 1000; %in milliseconds - NOTE THAT THIS MAY BE USED FOR THE CONCENTRATION CALCULATION INSTEAD OF T1 MAPS 
% r1_Cagent = 0.0037; %Relaxometry of Contrast Agent used
r1_Cagent = 0.0039; %Relaxometry of Contrast Agent used jkc

if DCE_OPTS(10)==1
   init_pars_fit = [log(0.15) 3 ];% INIT pars kinetic fit TOFTS
elseif DCE_OPTS(10)==2
   init_pars_fit = [log(0.15) 3 log(0.01)];% INIT pars kinetic fit Extende TOFTS
end

DCE_OPTStxt.t1 =  {'allS_' 'selS_' 'allV_' 'selV_' 'bV_' 'TE1_'  'R2starcor_' 'WB_' 'popAIF_'  'Tfts_' 'inpHC_' 'subT1_'  'FIG_' 'smth_'};
DCE_OPTStxt.t2 =  {'allS_' 'selS_' 'allV_' 'selV_' 'bV_' 'TE2_'  'R2starcor_' 'ROI_' 'subAIF_'  'eTfts_' 'inpSVD_' 'T1HC' num2str(T1_fixed) '_'  'FIG_' 'smth_'};
DCE_OPTStxt.t3 =  {'allS_' 'selS_' 'allV_' 'selV_' 'bV_' 'TE2_'  'R2starcor_' 'Smartmap_' 'subAIF_'  'eTfts_' 'inpSVD_' 'DCEbT1_'  'FIG_' 'smth_'};

% check to see if dicom parameters were read properly:
dicomissuecatchID=[]; 

% CREATE OUTPUT FOLDERS TREE:
% Create RUN storage folder with chosen DCE_OPTS in name:
txt1 =  [DCE_OPTStxt.t1{DCE_OPTS==1}];txt2 =  [DCE_OPTStxt.t2{DCE_OPTS==2}];txt3 =  [DCE_OPTStxt.t3{DCE_OPTS==3}];
RUN_dirname = [outputVisitDir 'RUN_' txt1 txt2 txt3]; 
RUN_dirname = regexprep(RUN_dirname,'__','_');
if RUN_dirname(end) == '_'
   RUN_dirname = RUN_dirname(1:end-1);
end

% Create MAT files storage folder within RUN folder:
MATfiles_dirname = [RUN_dirname '/MAT_FILES/']; 
if ~exist(MATfiles_dirname)
    mkdir(MATfiles_dirname);
end

% Create DCE_PROC MAPS files storage folder within RUN folder:
DCE_PROC_MAPSfiles_dirname = [RUN_dirname '/DCE_PROC/MAPS/']; 
if ~exist(DCE_PROC_MAPSfiles_dirname)
    mkdir(DCE_PROC_MAPSfiles_dirname);
end
% Create DCE_PROC MAPS/JPG files storage folder within MAPS folder:
DCE_PROC_MAPSJPGfiles_dirname = [RUN_dirname '/DCE_PROC/MAPS/JPG/']; 
if ~exist(DCE_PROC_MAPSJPGfiles_dirname)
    mkdir(DCE_PROC_MAPSJPGfiles_dirname);
end
% Create DCE_PROC MAPS/NII files storage folder within MAPS folder:
DCE_PROC_MAPSNIIfiles_dirname = [RUN_dirname '/DCE_PROC/MAPS/NII/']; 
if ~exist(DCE_PROC_MAPSNIIfiles_dirname)
    mkdir(DCE_PROC_MAPSNIIfiles_dirname);
end
% Create DCE_PROC STATS files storage folder within RUN folder:
DCE_PROC_STATSfiles_dirname = [RUN_dirname '/DCE_PROC/STATS/']; 
if ~exist(DCE_PROC_STATSfiles_dirname)
    mkdir(DCE_PROC_STATSfiles_dirname);
end
% Create TMP folder for leftover tests and anything else you want to save:
if ~exist([DCE_PROC_MAPSfiles_dirname(1:end-5) 'ALL_VOX/'])
    mkdir([DCE_PROC_MAPSfiles_dirname(1:end-5) 'ALL_VOX/']);
end


% END OF INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START DCE PROCESSING                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp ('=========================================================')

patID = 'DRO';


%%%% STEP 1 SET SCAN PARAMETERS
        % A 5 minute study is simulated, with injection of contrast agent occurring at 60 seconds.  
        % Thus, the total duration of imaging is 360 seconds.

        TR  = 5; % Repetition time in millisec! 
        alpha=30;alpha_rad = (pi/180).*alpha;
%         FR = 6; % SARA: seconds, but used in minutes in all scripts....
        FR = 0.5; %jkc

        FR_mins = FR/60;    % in minutes 
%         Total_scan_time_mins = 5;
        Total_scan_time_mins = 11; %jkc
        nr_of_frames =  Total_scan_time_mins/FR_mins+1; 

%         firstbaseline = 2; 
%         lastbaseline = 1/FR_mins; 
%         lpbs=lastbaseline+1;     %dynamic #
%% jkc

          firstbaseline = 2; 
          lastbaseline = 120; 
          lpbs=lastbaseline+1;     %dynamic #


        
%%%% STEP 2 % CALCULATE CONCENTRATION OVER TIME
        Step2_DRO_calculate_concentration(PreprocDataDir,patID,MATfiles_dirname,[],[],TR,T1_fixed,alpha_rad,r1_Cagent,firstbaseline,lastbaseline,DCE_OPTS);

        

%%%% STEP 3 % CALCULATE AIF   (depending on chosen option)     
        AIF = Step3_DRO_generateAIF(nr_of_frames,FR_mins,lpbs);
        
        
%%%% STEP 4 % Let the DCE kinetic modeling begin
        Step4_DRO_evaluateToftsModel(patID,MATfiles_dirname,PreprocDataDir,[],DCE_PROC_MAPSfiles_dirname,DCE_PROC_STATSfiles_dirname,DCE_OPTS,FR_mins,lpbs,AIF,init_pars_fit,firstbaseline,lastbaseline);
          

disp('Done! Have a great day.')

%% LOOK AT DATA:
cd /autofs/cluster/qtim/users/sara/DCE_DRO/RUN_selS_selV_TE1_R2starcor_popAIF_Tfts_inpHC_subT1_FIG_ROI/DCE_PROC/MAPS
load('DRO_KtransVeK2T1_MAPS.mat');

%%
figure,
imshow(Ktransmap,[0 1]);colormap jet;colorbar

figure,
imshow(Ve,[0 1]);colormap jet;colorbar

figure,
imshow(K2,[0 0.03]);colormap jet;colorbar


%% WHAT IT SHOULD BE
% x     y   Ktrans	Ve
% 0     10	0.01	0.01
% 0     20	0.02	0.01
% 0     30	0.05	0.01
% 0     40	0.1     0.01
% 0     50	0.2     0.01
% 0     60	0.35	0.01
% 10	10	0.01	0.05
% 10	20	0.02	0.05
% 10	30	0.05	0.05
% 10	40	0.1     0.05
% 10	50	0.2     0.05
% 10	60	0.35	0.05
% 20	10	0.01	0.1
% 20	20	0.02	0.1
% 20	30	0.05	0.1
% 20	40	0.1 	0.1
% 20	50	0.2 	0.1
% 20	60	0.35	0.1
% 30	10	0.01	0.2
% 30	20	0.02	0.2
% 30	30	0.05	0.2
% 30	40	0.1     0.2
% 30	50	0.2     0.2
% 30	60	0.35	0.2
% 40	10	0.01	0.5
% 40	20	0.02	0.5
% 40	30	0.05	0.5
% 40	40	0.1     0.5
% 40	50	0.2     0.5
% 40	60	0.35	0.5

