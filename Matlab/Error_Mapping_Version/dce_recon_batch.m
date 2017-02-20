% parpool(10)
% 
addpath C:\Users\azb22\Documents\Scripting\DCE_Testing\Matlab_QTIM_Script\rdir\
addpath C:\Users\azb22\Documents\Scripting\DCE_Testing\Matlab_QTIM_Script\DCEEstimation\code\
addpath C:\Users\azb22\Documents\Scripting\DCE_Testing\Matlab_QTIM_Script\NIfTI_20140122\

output = 'C:\Users\azb22\Documents\Scripting\DCE_Testing\Matlab_QTIM_Script\_error_';
% output = 'C:\Users\azb22\Documents\Scripting\DCE_Testing\DCE_Challenge\QIBA_v9_Tofts_Siemens_Orig\ans\toftsV6_'

input_file = 'C:\Users\azb22\Documents\Scripting\PBR\DCE\PBR_01\VISIT_01\dce1_mc_ss.nii.gz';
output_folder = 'C:\Users\azb22\Documents\Scripting\PBR\DCE\PBR_01\VISIT_01\';
% dcerecon_popAIF_volumes(input_file, output_folder, '')
% dcerecon_popAIF_tofts('C:\Users\azb22\Documents\Scripting\DCE_Testing\Matlab_QTIM_Script\Test_data\tofts_v9_5SNR.nii', output, '','')
dcerecon_popAIF_tofts('C:\Users\azb22\Documents\Scripting\DCE_Testing\Matlab_QTIM_Script\Test_data\gradient_toftsv6.nii', output, '','')
% dcerecon_popAIF_tofts('C:\Users\azb22\Documents\Scripting\DCE_Testing\Matlab_QTIM_Script\Test_data\tofts_v9.nii', output, '')
% dcerecon_popAIF_tofts('C:\Users\azb22\Documents\Scripting\DCE_Testing\Matlab_QTIM_Script\Test_data\tofts_v6.nii.gz', output, '','')

% dcerecon_popAIF_tofts('C:\Users\azb22\Documents\Scripting\DCE_Testing\DCE_Challenge\All_DROs\QIBA_v06_Tofts_RevB_Siemens\DICOM\20120109_090000dynamics018a001.nii.gz', output, '','')

% for dcmfolder = dir('C:\Users\azb22\Documents\Scripting\DCE_Testing\DCE_Challenge\QIBA_v9_Tofts_Siemens_Orig\*')'
% %    dicm2nii(dcmfolder.name, dcmfolder.name, 1)
% if length(dcmfolder.name) > 4
% foldername = char(strcat('C:\Users\azb22\Documents\Scripting\DCE_Testing\DCE_Challenge\QIBA_v9_Tofts_Siemens_Orig\', dcmfolder.name, '\DICOM'))
%    dcmfilename = char(strcat(foldername, '\dynamic.nii.gz'))
% %    t1mapname = char(strcat('C:\Users\azb22\Documents\Scripting\DCE_Testing\DCE_Challenge\All_DROs\T1_Maps\',dcmfolder.name(11:end),'\DICOM\T1Map.nii'));
% %    if ~(exist(t1mapname , 'file') == 2)
% %    t1mapname = char(strcat('C:\Users\azb22\Documents\Scripting\DCE_Testing\DCE_Challenge\All_DROs\T1_Maps\',dcmfolder.name(12:end),'\DICOM\T1Map.nii'));       
% %    end
%    outputpath = char(strcat('C:\Users\azb22\Documents\Scripting\DCE_Testing\DCE_Challenge\QIBA_v9_Tofts_Siemens_Orig\ans\', dcmfolder.name, '_T1STATIC_'));
%    dcerecon_popAIF_tofts(dcmfilename, outputpath, '', '')
% end
% end