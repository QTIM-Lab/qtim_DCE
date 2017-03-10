addpath ./rdir/
addpath ./DCEEstimation/code/
addpath ./NIfTI_20140122/

%For CED/NHX data...
directories_echo2 = [{'/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_01/MAPS/dce_mc_st_eco2.nii'},...
   {'/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_01/**/dce_mc_st_eco2.nii'},...
   {'/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_02/**/dce_mc_st_eco2.nii'},...
    {'/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_02/MAPS/dce_mc_st_eco2.nii'}];

directories_echo1 = [{'/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_01/MAPS/dce_mc_st_eco1.nii'},...
   {'/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_01/**/dce_mc_st_eco1.nii'},...
   {'/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_02/**/dce_mc_st_eco1.nii'},...
    {'/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_02/MAPS/dce_mc_st_eco1.nii'}];

directories_corrected = [{'/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_01/MAPS/dce_mc_st_corrected.nii'},...
   {'/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_01/**/dce_mc_st_corrected.nii'},...
   {'/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_02/**/dce_mc_st_corrected.nii'},...
    {'/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_02/MAPS/dce_mc_st_corrected.nii'}];

% For TMZ data...
% directories_echo1 = [{'/qtim2/users/data/TMZ/ANALYSIS/DCE/TMZ_05/VISIT_01/dce1_mc_ss.nii.gz'}];
% directories_echo2 = [{'/qtim2/users/data/TMZ/ANALYSIS/DCE/TMZ_05/VISIT_01/dce2_mc_ss.nii.gz'}];
% directories_corrected = [{'/qtim2/users/data/TMZ/ANALYSIS/DCE/TMZ_05/VISIT_01/dce_rstar2_corrected.nii.gz'}];

% This is a strange way to do this.. Not sure how to iterate through
% variables in Matlab.
directories = [directories_echo1, directories_echo2, directories_corrected];
directory_lengths = [0, length(directories_echo1), length(directories_echo2) + length(directories_echo2), length(directories_echo2) + length(directories_echo2) + length(directories_corrected)];
output_folder = [{'Echo1'},{'Echo2'},{'Corrected'}];

% Changes these parameters..
Destination_Path = '/home/abeers/Data/DCE_Package/Test_Results/';
Test_Code = '_Function_Test_';
Use_ROI = 1;
Use_T1Map = [0, 1];
Use_IndividualAIF = [0, 1];

gaussian_kernel_blur = [0,.2,.4,.6,.8,1,1.2];
PCA_levels = [0,1,2,3,4,5];
gd_dose = 1;
flip_angle = 10;
T1_tissue = 1500;
T1_blood = 1440;
TR = 6.8;
total_scan_time_seconds = 360;
bolus_arrival_time_seconds = 160;
hematocrit = .45;
relaxivity = .0045;
noise_threshold = [-1, 0, .01, .1];
T1_blur = gaussian_kernel_blur;
integration_method = {'recursive', 'conv'};
fitting_method = {'simplex', 'l-m'};
PCA_output = 0;
processes = 28;


for dce_type_num = 1
    for subdirectory = (directory_lengths(dce_type_num)+1):directory_lengths(dce_type_num+1)
        dce_files = rdir(char(directories(subdirectory)));
        for dce = dce_files'
            
            % Pre-Processing Options
            for i_gaussian_kernel_blur = gaussian_kernel_blur
            for i_PCA_levels = PCA_levels
            for i_noise_threshold = noise_threshold
            
            % Injection Parameter Options
            for i_gd_dose = gd_dose
            for i_total_scan_time_seconds = total_scan_time_seconds
            for i_flip_angle = flip_angle
            for i_T1_tissue = T1_tissue
            for i_T1_blood = T1_blood
            for i_TR = TR
            for i_bolus_arrival_time_seconds = bolus_arrival_time_seconds
            for i_hematocrit = hematocrit
            for i_relaxivity = relaxivity
                
            % Fitting Parameter Options
            for i_integration_method = integration_method
            for i_fitting_method = fitting_method
            
            % Input Parameter Map Options
            for i_Use_T1Map = Use_T1Map
            for i_Use_IndividualAIF = Use_IndividualAIF
            
            % Make Output Folder
            if ~exist(strcat(Destination_Path, char(output_folder(dce_type_num))), 'dir')
                mkdir(strcat(Destination_Path, char(output_folder(dce_type_num))));
            end
            
            % Print Current File
            disp(dce.name)
            
            % Get Patient/Visit Prefix, Find Output Path
            [pathstr, name, ext] = fileparts(dce.name);
            temp = strsplit(dce.name, filesep);
            outputpath = strcat(temp(end-3), '_', temp(end-2));
            outputpath = strcat(Destination_Path, char(output_folder(dce_type_num)), filesep, char(outputpath));
            
            % Get Parameter Map Directories
            autoAIF_dir = char(strjoin(temp(1:end-1),filesep));
            ROI_dir = char(strjoin(temp(1:end-2),filesep));
            ROI_dir = char(strcat(ROI_dir, '/ROISTATS/T1AxialPost/'));

            try
                
                % Individual AIF Check
                if (i_Use_IndividualAIF == 1)
                found_AIF = rdir(char(strcat(autoAIF_dir, '/*NORDIC*')));
                autoAIF_file = found_AIF.name;
                outputpath = char(strcat(outputpath, '_autoAIF_'));
                else
                    autoAIF_file = '';
                end

                % T1Map Check
                if (i_Use_T1Map == 1)
                found_T1Map = rdir(char(strcat(autoAIF_dir, '/*T1in*')));
                T1Map_file = found_T1Map.name;
                outputpath = char(strcat(outputpath, '_t1map_'));
                else
                    T1Map_file = '';
                end

                % ROI Check
                if (Use_ROI == 1)
                found_ROI = char(strcat(ROI_dir, 'rT1AxialPostROI.nii'));
                ROI_file = found_ROI;
                else
                    ROI_file = '';
                end

                % Integration Method Check
                if length(integration_method) > 1
                outputpath = char(strcat(outputpath,'_',i_integration_method{1},'_'));
                end
                % Fitting Method Check
                if length(fitting_method) > 1
                outputpath = char(strcat(outputpath,'_',i_fitting_method{1},'_'));
                end
                % Blur Check
                if length(gaussian_kernel_blur) > 1
                outputpath = char(strcat(outputpath,'_blur_',num2str(i_gaussian_kernel_blur),'_'));
                end                
                % PCA Check
                if length(PCA_levels) > 1
                outputpath = char(strcat(outputpath,'_pca_',num2str(i_PCA_levels),'_'));
                end                      
                % Threshold Check
                if length(noise_threshold) > 1
                outputpath = char(strcat(outputpath,'_threshold_',num2str(i_noise_threshold),'_'));
                end     
                
                % outputpath = char(strcat(outputpath,'_kernel_',num2str(kernel_blur),'_'));
                outputpath = char(strcat(outputpath, Test_Code));
                
                % Check if this file has already been run.
                if exist(strcat(outputpath, 'ktrans.nii.gz'))
                   continue;
                end

                DCE_Processing(...
                dce.name,... %input_file
                outputpath,... %output_path
                ROI_file, ... %mask_file
                autoAIF_file,... %provided_AIF
                T1Map_file,... %T1_map
                i_gaussian_kernel_blur,... %gaussian_kernel_blur
                i_PCA_levels,... %PCA_levels
                i_gd_dose,... %gd_dose
                i_flip_angle,... %flip_angle
                i_T1_tissue,... %T1_tissue
                i_T1_blood,... %T1_blood
                i_TR,... %TR
                i_total_scan_time_seconds,... % total_scan_time_seconds
                i_bolus_arrival_time_seconds,... % bolus_arrival_time_seconds
                i_hematocrit,...% hematocrit
                i_relaxivity,... % relaxivity
                i_noise_threshold,... % noise_threhold
                T1_blur,... % T1_blur
                i_integration_method,... % integration method
                i_fitting_method,... % fitting method
                PCA_output,... % PCA_output
                processes) % processes
        
            catch
                disp('Error! Cancelling Processing');
            end
            
            end
            end
            end
            end
            end
            end
            end
            end
            end
            end 
            end
            end
            end
            end
            end
            end
        end
    end
end