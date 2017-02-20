%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   STEP 4a                       %
%       Preparation for the kinetic modeling      %
%     of the DCE MRI data and looping over all    %
%               masked voxels                     %
%                                                 %
%                                                 %
% Sara Maria Sprinkhuizen - March 2013         	  %
% based on code of                                %
% Kim Mouridsen  <kimm@nmr.mgh.harvard.edu>       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
function [Ktrans, Ve, K2, Vp, est_conc_mat, obs_conc_mat, Curves_N_Pars_allmaskvoxels] = Step4a_DRO_Kinmap(conc_4D,FR_mins,lpbs,AIF,init_pars_fit,mask,MDL,DCE_OPTS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT  
% - conc_4D        : 4D contrast agent concentration volume
% - FR_mins        : Frame rate/dynamic scan time in minutes
% - lpbs           : Frame number of frame @ which CA bolus arrives
% - AIF            :  Arterial Input Function
% - init_pars_fit  : Initial fitting parameters for Ktrans, Ve, (Vp) 
% - mask           : mask within which kinetic modeling is done
% - MDL            : 1= Tofts, 2 = Extended Tofts 
% - DCE_OPTS       : user chosen processing options
%
%
% OUTPUT
% - Ktrans
% - Ve
% - K2
% - Vp
% - est_conc_mat
% - obs_conc_mat
% - Curves_N_Pars_allmaskvoxels : contains the observed and estimated (fitted) curves of ALL masked voxels
%                                 as well as all individual Ktrans, Ve, Kep values per masked voxel
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
nrow = size(conc_4D,1); ncol = size(conc_4D,2); 
ind = find(conc_4D(:,:,1)>-10);
[Ix, Iy] = ind2sub([nrow ncol], ind);

conc_4D_nobaseline = squeeze(conc_4D(:,:,(lpbs+1):end));
AIF_nobaseline = AIF((lpbs+1):end);

%% Ktrans Ve Kep (Vp) using (extended) Toft's model
Ktrans = zeros(nrow,ncol);
K2 = zeros(nrow,ncol);
Ve = zeros(nrow,ncol);
Vp = zeros(nrow,ncol);

est_conc_mat = zeros(length(Ix),length(AIF_nobaseline));
obs_conc_mat = zeros(length(Ix),length(AIF_nobaseline));

% SARA ADDED AS TEST:
Curves_N_Pars_allmaskvoxels = [];matcount = 1;
obs_conc_alldyns = zeros(length(Ix),length(AIF));

disp('   Starting the DCE kinetic modeling')
for i=1:length(Ix)
  if mod(i,100)==0
    disp(['   Voxel ' num2str(i) ' / ' num2str(length(Ix))]);
  end
  % Get the CA concentration curve in this voxel
  concvec = squeeze(conc_4D_nobaseline(Ix(i),Iy(i),:));
 
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
if MDL == 1  % TOFTS

%%%%%%%%%% FITTING BEGINS
  [kin_fit,est_conc] = Step4b_Simplex_Fit_Kinetic_Tofts(concvec,AIF_nobaseline',init_pars_fit,FR_mins);
  
%%%%%%%%%% Processing of outcome of fitting
  Ktrans(Ix(i), Iy(i)) = (exp(kin_fit(1)));
  Ve(Ix(i), Iy(i)) = 1/(1+(exp(-kin_fit(2))));
  K2(Ix(i), Iy(i)) = Ktrans(Ix(i), Iy(i)) /  Ve(Ix(i), Iy(i));
  est_conc_mat(i,:) = est_conc;
  obs_conc_mat(i,:) = concvec;

  % Sara added as test output:
  obs_conc_alldyns = squeeze(conc_4D(Ix(i),Iy(i),:));
  est_conc_alldyns = [zeros(1,lpbs) est_conc]; 
  Curves_N_Pars_allmaskvoxels(matcount).est_conc = est_conc_alldyns;
  Curves_N_Pars_allmaskvoxels(matcount).obs_conc = obs_conc_alldyns;
  Curves_N_Pars_allmaskvoxels(matcount).Ktrans =  (exp(kin_fit(1)));
  Curves_N_Pars_allmaskvoxels(matcount).Ve = 1/(1+(exp(-kin_fit(2))));
  Curves_N_Pars_allmaskvoxels(matcount).Kep = (exp(kin_fit(1))) / (1/(1+(exp(-kin_fit(2)))));
  matcount = matcount+1;
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
elseif MDL == 2  % Extended TOFTS
%%%%%%%%%% FITTING BEGINS
  [kin_fit,est_conc] = Step4b_Simplex_Fit_Kinetic_Extended_Tofts(concvec,AIF_nobaseline',init_pars_fit,FR_mins);
  
%%%%%%%%%% Processing of outcome of fitting
  Ktrans(Ix(i), Iy(i)) = (exp(kin_fit(1)));
  Ve(Ix(i), Iy(i)) = 1/(1+(exp(-kin_fit(2))));
  Vp(Ix(i), Iy(i)) = exp(kin_fit(3));
  K2(Ix(i), Iy(i)) = Ktrans(Ix(i), Iy(i)) /  Ve(Ix(i), Iy(i));
  
  est_conc_mat(i,:) = est_conc;
  obs_conc_mat(i,:) = concvec;

  % Sara added as test output:
  obs_conc_alldyns = squeeze(conc_4D(Ix(i),Iy(i),:));
  est_conc_alldyns = [zeros(1,lpbs) est_conc]; 
  Curves_N_Pars_allmaskvoxels(matcount).est_conc = est_conc_alldyns;
  Curves_N_Pars_allmaskvoxels(matcount).obs_conc = obs_conc_alldyns;
  Curves_N_Pars_allmaskvoxels(matcount).Ktrans =  (exp(kin_fit(1)));
  Curves_N_Pars_allmaskvoxels(matcount).Ve = 1/(1+(exp(-kin_fit(2))));
  Curves_N_Pars_allmaskvoxels(matcount).Kep = (exp(kin_fit(1))) / (1/(1+(exp(-kin_fit(2)))));
  Curves_N_Pars_allmaskvoxels(matcount).Vp = exp(kin_fit(3));
  matcount = matcount+1;
    

end

end


