function estconcen4D = estiamteConcentrationsfromObservedConcentrations(obsConcenImg4DName, estConcenImg4DName);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program reads the observed concentration image (4D) and outputs the estimated concentration image (4D).
%
% Usage:
% 
% estconcen4D = estimateConcentrationsfromObservedConcentrations(obsConcenImg4DName, estConcenImg4DName)
%
% where
%
%          obsConcenImg4DName --- [input]  the name of the obsered concentration file (a 4D image)
%          estConcenImg4DName --- [output] the name of the estimated concentration file (a 4D image)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (nargin<2)
	fprintf('Error: Not enought input arguments!\n');
	fprintf('Please specify an input 4D (observed) concentration file name and an output 4D (estimated) concentration file name!\n\n');
	fprintf('Or, run "help estimateConcentrationsfromObservedConcentrations" to see how to use.\n\n');
	return;
end


% library
addpath /autofs/cluster/qtim/users/you2/DCEEstimation/code              % dce package
addpath /autofs/cluster/mi2b2/4/you2/Utilities/NIfTI_matlab_20140122/   % nifti read

% begin prompt
fprintf('This program reads the observed concentration image and outputs the estimated concentration\n');


% read image and parse its dimension
concen4D=load_untouch_nii(obsConcenImg4DName);
fprintf('The observed concentration image has been read: %s\n', obsConcenImg4DName);
xSize = size(concen4D.img,1);
ySize = size(concen4D.img,2);
zSize = size(concen4D.img,3);
tSize = size(concen4D.img,4);
uSize = size(concen4D.img,5);


% number of effective frames and where it starts
nbaseline=11;                   % the number of baseline frames
nSize = max(tSize, uSize);      % number of frames
neSize= nSize-nbaseline+1;      % number of effective frames (excluding all but the last baseline frame)
tStart= min(tSize, nbaseline);
tEnd  = tSize;
uStart= min(uSize, nbaseline);
uEnd  = uSize;

% QC 
if (nSize<=1)
    fprintf('Error: the input concentration image is not 4D! Program exits. \n');
    return;
end



% AIF and init_parameters are the same everywhere in the image
FR=6;                          % seconds, time between frames, in seconds
FR_mins=FR/60;                 % FR in minutes
lpbs=1;                        % last point for baseline, set to 1 because we are excluding all but the last baseline frames
init_parameters=[log(0.15) 3]; % initial guess for Ktrans and Ve
AIF=generateAIF(neSize, FR, lpbs);  % generate an AIF of the effective number of frames



% prepare for the output concentration image
estconcen4D=concen4D;
estconcen4D.hdr.dime.dim(5)=nSize;
estconcen4D.hdr.dime.dim(6)=1;
estconcen4D.img=zeros(xSize,ySize,zSize,nSize);

% fit the concentration for each voxel in the image
for k=1:zSize
  fprintf('slice %d...\n', k);
  for i=1:xSize
    for j=1:ySize
	%fprintf('fitting the concentration for voxel (%d, %d, %d)\n', i, j, k);
	[kin_para est_concen] = Step4b_Simplex_Fit_Kinetic_Tofts(squeeze(concen4D.img(i,j,k,tStart:tEnd,uStart:uEnd)), AIF(nbaseline:end), init_parameters, FR_mins);
	estconcen4D.img(i,j,k,nbaseline:end) = est_concen';
    end
  end
end

% save the 4D estimated concentration image as the output
save estconcen4D estconcen4D;   % save .mat to backup (this line should be removed once the saving into image line below functions well)
save_untouch_nii(estconcen4D, estConcenImg4DName);
fprintf('the estimated concentration image has been saved into %s\n\n', estConcenImg4DName);
fprintf('\nthe program exits successfully.\n\n');

