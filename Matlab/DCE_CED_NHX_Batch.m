% Only run this line once..
%parpool(28)

addpath /home/abeers/Projects/Matlab_DCE/Matlab_QTIM_Script/rdir/
addpath /home/abeers/Projects/Matlab_DCE/Matlab_QTIM_Script/DCEEstimation/code/
addpath /home/abeers/Projects/Matlab_DCE/Matlab_QTIM_Script/NIfTI_20140122/

%For CED/NHX data...
% directories_echo2 = [{'/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_01/MAPS/dce_mc_st_eco2.nii'},...
%    {'/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_01/**/dce_mc_st_eco2.nii'},...
%    {'/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_02/**/dce_mc_st_eco2.nii'},...
%     {'/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_02/MAPS/dce_mc_st_eco2.nii'}];
% 
% directories_echo1 = [{'/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_01/MAPS/dce_mc_st_eco1.nii'},...
%    {'/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_01/**/dce_mc_st_eco1.nii'},...
%    {'/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_02/**/dce_mc_st_eco1.nii'},...
%     {'/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_02/MAPS/dce_mc_st_eco1.nii'}];
% 
% directories_corrected = [{'/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_01/MAPS/dce_mc_st_corrected.nii'},...
%    {'/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_01/**/dce_mc_st_corrected.nii'},...
%    {'/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_02/**/dce_mc_st_corrected.nii'},...
%     {'/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_02/MAPS/dce_mc_st_corrected.nii'}];

% For TMZ data...
directories_echo1 = [{'/qtim2/users/data/TMZ/ANALYSIS/DCE/TMZ_05/VISIT_01/dce1_mc_ss.nii.gz'}];
directories_echo2 = [{'/qtim2/users/data/TMZ/ANALYSIS/DCE/TMZ_05/VISIT_01/dce2_mc_ss.nii.gz'}];
directories_corrected = [{'/qtim2/users/data/TMZ/ANALYSIS/DCE/TMZ_05/VISIT_01/dce_rstar2_corrected.nii.gz'}];

directories = [directories_echo1, directories_echo2, directories_corrected];
directory_lengths = [0, length(directories_echo1), length(directories_echo2) + length(directories_echo2), length(directories_echo2) + length(directories_echo2) + length(directories_corrected)];
output_folder = [{'Echo1'},{'Echo2'},{'Corrected'}];

% Changes these parameters..
roi = 0;
Flip_Angle = 10;
Static_T1 = 1000;
Destination_Path = '/home/abeers/Projects/Matlab_DCE/Matlab_QTIM_Script/Test_Results/';
Blur = [.65];
PCA_Components = [0];
Dose = .5;
Threshold = .01;
    
for dce_type_num = 1
    for subdirectory = (directory_lengths(dce_type_num)+1):directory_lengths(dce_type_num+1)
        directories(subdirectory)
        dce_files = rdir(char(directories(subdirectory)));
        for dce = dce_files'
            %if strfind(char(dce.name),'OLD') == 0
            %   continue 
            %end
            for kernel_blur = Blur
                for PCA = PCA_Components
                    for t1_status = 1
                        for aif_status = 1
                    
                    if ~exist(strcat(Destination_Path, char(output_folder(dce_type_num))), 'dir')
                        mkdir(strcat(Destination_Path, char(output_folder(dce_type_num))));
                    end
                   
                        
                            dce.name
                            temp = strsplit(dce.name, '/');
                            paths = length(temp);
                            outputpath = strcat(temp(paths-2), '_', temp(paths-1));
%                             if (roi == 1)
%                             outputpath = strcat('/qtim2/users/data/NHX/ANALYSIS/DCE/DCE_KTRANS_ANALYSIS/ROI_ONLY/', char(output_folder(dce_type_num)), '/', char(outputpath),'_');
%                             else
%                             outputpath = strcat('/qtim2/users/data/NHX/ANALYSIS/DCE/DCE_KTRANS_ANALYSIS/FULL_VOL/', char(output_folder(dce_type_num)), '/', char(outputpath),'_');
%                             end
%                             outputpath = strcat('/home/abeers/Requests/Jayashree/TMZ_Data/TMZ_DCE_KTRANS_Andrew_Code/', char(output_folder(dce_type_num)), '/', char(outputpath),'_');
                            outputpath = strcat(Destination_Path, char(output_folder(dce_type_num)), '/', char(outputpath),'_');
                            autoAIF_dir = char(strjoin(temp(1:end-1),'/'));
                            
                            ROI_dir = char(strjoin(temp(1:end-2),'/'));
                            ROI_dir = char(strcat(ROI_dir, '/ROISTATS/T1AxialPost/'));

                            suffix = '';

                            %try
                                if (aif_status == 2)
                                found_AIF = rdir(char(strcat(autoAIF_dir, '/*NORDIC*')))
                                autoAIF_file = found_AIF.name
                                outputpath = char(strcat(outputpath, '_autoAIF_'))
                                else
                                    autoAIF_file = '';
                                end

                                if (t1_status == 2)
                                found_T1Map = rdir(char(strcat(autoAIF_dir, '/*T1in*')))
                                T1Map_file = found_T1Map.name
                                outputpath = char(strcat(outputpath, '_t1map_'))
                                else
                                    T1Map_file = '';
                                end

                                if (roi == 1)
                                found_ROI = char(strcat(ROI_dir, 'rT1AxialPostROI.nii'))
                                ROI_file = found_ROI
                                else
                                    ROI_file = '';
                                end

                                outputpath = char(strcat(outputpath,'_kernel_',num2str(kernel_blur),'_'));
                                outputpath = char(strcat(outputpath,'_PCA_',num2str(PCA),'_'));
                                outputpath = char(strcat(outputpath,'_dose_',num2str(Dose),'_'));
                                outputpath = char(strcat(outputpath,'_threshold_',num2str(Threshold),'_'));
                                
                                dce.name;
                                ROI_file;
                                T1Map_file;
                                autoAIF_file;
                                t1_status;
                                aif_status;
                                roi;

                                dcerecon_popAIF_tofts(dce.name, outputpath, ROI_file, autoAIF_file, T1Map_file, kernel_blur, PCA, Dose, Flip_Angle, Static_T1, Threshold);
                            %catch
                            %    a = 0;
                            %end
                        end
                    end 
                end
            end
        end
    end
end