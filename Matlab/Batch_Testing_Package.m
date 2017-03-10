% Only run this line once..
%parpool(28)

addpath /home/abeers/Projects/DCE_Package/Matlab/rdir/
addpath /home/abeers/Projects/DCE_Package/Matlab/DCEEstimation/code/
addpath /home/abeers/Projects/DCE_Package/Matlab/NIfTI_20140122/

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

directories = [directories_echo1, directories_echo2, directories_corrected];
directory_lengths = [0, length(directories_echo1), length(directories_echo2) + length(directories_echo2), length(directories_echo2) + length(directories_echo2) + length(directories_corrected)];
output_folder = [{'Echo1'},{'Echo2'},{'Corrected'}];

% Changes these parameters..
roi = 1;
Flip_Angle = 10;
Static_T1 = 1000;
Destination_Path = '/home/abeers/Projects/DCE_Package/Matlab/Test_Results/';
Blur = [.65];
PCA_Components = [0];
Dose = .5;
Threshold = .01;
Test_Code = '_Bad_Corrected_HighT1'

for dce_type_num = 3
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
                            outputpath = strcat(temp(paths-3), '_', temp(paths-2));
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

                            try
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

                                % outputpath = char(strcat(outputpath,'_kernel_',num2str(kernel_blur),'_'));
                                % outputpath = char(strcat(outputpath,'_PCA_',num2str(PCA),'_'));
                                % outputpath = char(strcat(outputpath,'_dose_',num2str(Dose),'_'));
                                % outputpath = char(strcat(outputpath,'_threshold_',num2str(Threshold),'_'));
                                outputpath = char(strcat(outputpath, Test_Code));

                                if exist(strcat(outputpath, 'ktrans.nii.gz'))
                                   continue;
                                end

                                dce.name;
                                ROI_file;
                                T1Map_file;
                                autoAIF_file;
                                t1_status;
                                aif_status;
                                roi;

                                % dcerecon_popAIF_tofts(dce.name, outputpath, ROI_file, autoAIF_file, T1Map_file, kernel_blur, PCA, Dose, Flip_Angle, Static_T1, Threshold);
                                DCE_Processing(...
                                dce.name,... %input_file
                                outputpath,... %output_path
                                ROI_file, ... %mask_file
                                '',... %provided_AIF
                                '',... %T1_map
                                0,... %gaussian_kernel_blur
                                2,... %PCA_levels
                                1,... %gd_dose
                                10,... %flip_angle
                                1500,... %T1_tissue
                                1440,... %T1_blood
                                6.8,... %TR
                                360,... % total_scan_time_seconds
                                160,... % bolus_arrival_time_seconds
                                .45,...% hematocrit
                                .0045,... % relaxivity
                                .01,... % noise_threhold
                                0,... % T1_blur
                                1,... % parallel_status
                                0,... % PCA_output
                                28) % processes
                            catch
                               a = 0;
                            end
                        end
                    end 
                end
            end
        end
    end
end