% parpool(28)

addpath /home/abeers/Projects/Matlab_DCE/Matlab_QTIM_Script/rdir/
addpath /home/abeers/Projects/Matlab_DCE/Matlab_QTIM_Script/DCEEstimation/code/
addpath /home/abeers/Projects/Matlab_DCE/Matlab_QTIM_Script/NIfTI_20140122/

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

directories_echo2 = [{'/qtim/users/data/TMZ/ANALYSIS/DCE/**/**/dce_mc_st_eco2.nii'}];
directories_echo1 = [{'/qtim/users/data/TMZ/ANALYSIS/DCE/**/**/dce_mc_st_eco2.nii'}];
directories_corrected = [{'/qtim/users/data/TMZ/ANALYSIS/DCE/**/**/dce_mc_st_eco2.nii'}];

directories = [directories_echo2, directories_echo1, directories_corrected];
output_folder = [{'Echo1'},{'Echo2'},{'Corrected'}];

directory = '/home/abeers/Projects/DCE_Challenge/NIIs/08*DCE-2/*.nii';

dce_files = rdir(char(directory));


aif_status=2;
t1_status=2;

roi=1;

% for dce = dce_files'
% %     try
%     dce.name
%     temp = strsplit(dce.name, '/');
%     paths = length(temp);
%     outputpath = temp(paths-1);
%     outputpath = strcat('/home/abeers/Projects/DCE_Challenge/Results/', char(outputpath),'_');
%     
%     if aif_status==2
%        AIF_file = char('/home/abeers/Projects/DCE_Challenge/Patient_Population_AIF.txt');
%     end
%     
%     if t1_status==2
%         temp = strsplit(char(temp(paths-1)), '_');
%         T1_path = char(strcat(temp(2), '_', temp(3)));
%         if strcmp(T1_path(1),'0')
%             T1_path = T1_path(2:end)
%         end
%         T1_file = char(strcat('/home/abeers/Projects/DCE_Challenge/Angle_Niftis/', T1_path, '_T1Map.nii.gz'));
%     end
%     
%     AIF_file
%     T1_file
%     outputpath
%     
%     dcerecon_popAIF_tofts(dce.name, outputpath, '', AIF_file, T1_file, .65, 0, 1)
% %     catch
% %         a=0;
% %     end
%     
% end


% tmz_files = rdir('/qtim2/users/xda/Test_QTIM_Pipeline/TMZ_DATA/ANALYSIS/DCE/**/**/dce1_mc_ss.nii.gz');

% for i = 0:.1:3
%     AIF=generateAIF(250*2,6.86,61);
%     postInjection= AIF(61:(500-61));
%     newAIF = AIF;
%     for p = 61:2:(500-61)
%         newAIF(61 + (p-61)/2) = (AIF(p) + AIF(p+1))/2;
%     end
%     
%     newAIF = newAIF(1:250);
%     newAIF = newAIF./i;
%    
%     
%     oldAIF = generateAIF(250, 6.86, 61);
%    
%     plot(newAIF(61:65))
%     
%     tester = 2*trapz(newAIF(61:70)) - trapz(oldAIF(61:70))
%     i
%     
%     end
    
% DCE_Challenge_Patients = rdir('/home/abeers/Projects/DCE_Challenge/**/*DCE*.nii')
% 
% for visit = DCE_Challenge_Patients'
%     
%     splitfile = strsplit(visit,'/');
%     ID = char(splitfile(end-1));
%        
%     if strcmp(ID(1), 0) == 1
%         match_string = ID(5:7);
%     else
%         match_string = ID(4:7);
%     end
%     
%     T1MAP = strcat('/home/abeers/Projects/DCE_Challenge/Angle_Niftis/', match_string, '_T1Map.nii.gz')
%     
%     outputpath = strcat('/home/abeers/Projects/DCE_Challenge/Results/', match_string, '_MGH_')
%     
%     dcerecon_popAIF_tofts(inputpath, outputpath, '', '', T1Map, .65, 0, 1)
% 
% end

% trapz(obs_conc)/trapz(gd_AIF);
% 
% for patient = 1:8
%     for visit = 1:3
%     
%     inputpath = strcat('/qtim2/users/xda/Test_QTIM_Pipeline/TMZ_DATA/ANALYSIS/DCE/TMZ_0', num2str(patient), '/VISIT_0', num2str(visit), '/dce1_mc_ss.nii.gz');
%     inputpath = char(inputpath)
%  
%     
%     if (p > 2)
%         dose = .5;
%     else
%         dose = 1;
%         end
%     
%     outputpath = strcat('/home/abeers/Requests/Jayashree/TMZ_KTrans_Dose_Adjusted/TMZ_0', num2str(patient), '-VISIT_0', num2str(visit),'_registered_');
%     outputpath = char(outputpath)
%         
%     if exist(inputpath, 'file') == 2
% %     dcerecon_popAIF_tofts(inputpath, outputpath, '', '', '', .65, 0, dose)
%     end
%         
%     end 
% end

% for kernel_blur = 0:1.35:1.35
%     for PCA = 
for kernel_blur = .65
    for PCA = 0
        for t1_status = 1:2
            for aif_status = 2
                for dce_type_num = 1:3
                    for subdirectory = 1:4
                        if dce_type_num == 1
                        dce_files = rdir(char(directories_echo1(subdirectory)));
                        elseif dce_type_num == 2
                        dce_files = rdir(char(directories_echo2(subdirectory)));    
                        else
                        dce_files = rdir(char(directories_corrected(subdirectory)))
                        end
                        
                        for dce = dce_files'
                            dce.name
                            temp = strsplit(dce.name, '/');
                            paths = length(temp);
                            outputpath = strcat(temp(paths-3), '_', temp(paths-2));
%                             if (roi == 1)
%                             outputpath = strcat('/qtim2/users/data/NHX/ANALYSIS/DCE/DCE_KTRANS_ANALYSIS/ROI_ONLY/', char(output_folder(dce_type_num)), '/', char(outputpath),'_');
%                             else
%                             outputpath = strcat('/qtim2/users/data/NHX/ANALYSIS/DCE/DCE_KTRANS_ANALYSIS/FULL_VOL/', char(output_folder(dce_type_num)), '/', char(outputpath),'_');
%                             end
                            outputpath = strcat('/home/abeers/Requests/Jayashree/TMZ_Data/TMZ_DCE_KTRANS_Andrew_Code/', char(output_folder(dce_type_num)), '/', char(outputpath),'_');
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

                                outputpath = char(strcat(outputpath,'_kernel_',num2str(kernel_blur),'_'));
                                outputpath = char(strcat(outputpath,'_PCA_',num2str(PCA),'_'))

                                dce.name
                                ROI_file
                                T1Map_file 
                                autoAIF_file
                                t1_status
                                aif_status
                                roi;

                                dcerecon_popAIF_tofts(dce.name, outputpath, ROI_file, autoAIF_file, T1Map_file, kernel_blur, PCA, 1)
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
% dce_files = rdir('/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_01/**/dce_mc_st_eco1.nii');
% 
%     for dce = dce_files'
%     dce.name;
%     temp = strsplit(dce.name, '/');
%     temp = strcat(temp(8), '_', temp(9));
%     outputpath = strcat('/qtim2/users/data/NHX/ANALYSIS/DCE/DCE_KTRANS_ANALYSIS/', char(temp),'_');
%     try
%     dcerecon_popAIF_tofts(dce.name, outputpath, '')
%     catch
%     a=0;
%     end
%     end 
% 
% dce_files = rdir('/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_01/MAPS/dce_mc_st_eco1.nii');
% 
%     for dce = dce_files'
%     dce.name;
%     temp = strsplit(string(dce.name), '/');
%     temp = strcat(temp(9), '_', temp(10),'_');
%     outputpath = char(strcat('/qtim2/users/data/NHX/ANALYSIS/DCE/DCE_KTRANS_ANALYSIS/',char(temp)))
%         try
%     dcerecon_popAIF_tofts(dce.name, outputpath, '')
%     catch
%     a=0;
%     end
%     end
%     
% dce_files = rdir('/qtim2/users/data/NHX/ANALYSIS/DCE/**/VISIT_02/**/dce_mc_st_eco1.nii');
% 
%     for dce = dce_files'
%     dce.name;
%     temp = strsplit(dce.name, '/');
%     temp = strcat(temp(8), '_', temp(9));
%     outputpath = strcat('/qtim2/users/data/NHX/ANALYSIS/DCE/DCE_KTRANS_ANALYSIS/', char(temp),'_')
%         try
%     dcerecon_popAIF_tofts(dce.name, outputpath, '')
%     catch
%     a=0;
%     end
%     end
%     
% dce_files = rdir('/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/**/VISIT_02/MAPS/dce_mc_st_eco1.nii');
% autoAIF = 'true';
% 
%     for dce = dce_files'
%     dce.name
%     temp = strsplit(dce.name, '/');
%     outputpath = strcat(temp(9), '_', temp(10));
%     outputpath = strcat('/qtim2/users/data/NHX/ANALYSIS/DCE/DCE_KTRANS_ANALYSIS/', char(outputpath),'_')
%     autoAIF_dir = char(strjoin(temp(1:end-1),'/'))
%     found_AIF = rdir(char(strcat(autoAIF_dir, '/*AIF*')))
%     autoAIF_file = found_AIF.name
%     
%     found_T1Map = rdir(char(strcat(autoAIF_dir, '/*T1in*')))
%     T1Map_file = found_T1Map.name
%   
%     dcerecon_popAIF_tofts(dce.name, outputpath, '', '', T1Map_file)
%     end    
    
% input_file = '/home/abeers/Junk/PBR_05_V02_mc_ss.nii.gz'
% input_file = '/home/abeers/Junk/P5V2_no_mc_ss.nii.gz'
% input_file = '/qtim2/users/data/PBR/ANALYSIS/DCE/PBR_05/VISIT_01/dce1_mc_ss.nii.gz';
% output_folder = '/home/abeers/Junk/';
% dcerecon_popAIF_tofts('C:\Users\azb22\Documents\Scripting\DCE_Testing\Matlab_QTIM_Script\Test_data\tofts_v9_5SNR.nii', output, '')
% dcerecon_popAIF_tofts('C:\Users\azb22\Documents\Scripting\DCE_Testing\Matlab_QTIM_Script\Test_data\tofts_v6.nii.gz', output_folder, '')
% dcerecon_popAIF_tofts('C:\Users\azb22\Documents\Scripting\DCE_Testing\Matlab_QTIM_Script\Test_data\tofts_v9.nii', output, '')
% dcerecon_popAIF_tofts(input_file, output_folder, '')