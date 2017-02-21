function DCE_Processing(...
input_file,...
output_path,...
mask_file, ...
provided_AIF,...
T1_map,...
gaussian_kernel_blur,...
PCA_levels,...
gd_dose,...
flip_angle,...
T1_tissue,...
T1_blood,...
TR,...
total_scan_time_seconds,...
bolus_arrival_time_seconds,...
hematocrit,...
relaxivity,...
noise_threshold,...
T1_blur,...
parallel_status,...
PCA_output,...
processes)

% Set Default Parameters
if (nargin==2)
	mask_file='';
    provided_AIF='';
elseif (nargin==3)
    mask_file='';
end

% Initialize useful function variables.
alpha=flip_angle;alpha_rad = (pi/180).*alpha;
first_baseline = 1;
R1_pre = 0;
T1_mapping = 0;
convert_AIF = 0;

% Load mask if available.
if ~(isempty(mask_file))
    if ~(exist(mask_file))
        error('Mask filename does not exist. Aborting processing.');
    else
        mask = load_untouch_nii(mask_file);
        mask_img = mask.img;
    end
end

% If input file exists, begin with processing.
if ~(exist(input_file))
    error('Input DCE Volume does not exist. Aborting processing.');
else
    dce4D=load_untouch_nii(input_file);
    xsize=size(dce4D.img,1);
    ysize=size(dce4D.img,2);
    zsize=size(dce4D.img,3);
    tsize=size(dce4D.img,4);
    
    % Check for T1 Map.
    % note: make this work irrespective of dimensions
    if ~(isempty(T1_map))
        T1_map_file = load_untouch_nii(T1_map);
        T1_map_img = T1_map_file.img;
        T1_mapping = 1;

        if (gaussian_kernel_blur > 0)
            for z = 1:size(T1_map_img,3)
                T1_map_img(:,:,z) = imgaussfilt(T1_map_img(:,:,z), gaussian_kernel_blur);
            end
        end

        R1_pre = 1./T1_map_img;
    else
        R1_pre = (1./T1_tissue).*ones([xsize ysize zsize]);
    end
    
    time_interval_seconds = total_scan_time_seconds/tsize;
    time_interval_mins = time_interval/60;
    last_baseline = round(bolus_arrival_time_seconds/time_interval_seconds);
    bolus_start = last_baseline + 1;

    % note: make this work irrespective of dimensions
    ktransmap=zeros(xsize,ysize,zsize,tsize);
    vemap=zeros(xsize,ysize,zsize,tsize);
    aucmap=zeros(xsize,ysize,zsize,tsize);

    % Parse provided AIF. This presumes the provided AIF is already
    % converted to signal.
    % note: perhaps make a flag to optionally set signal conversion on AIF.
    % note: be sure to provide specific instructions on how to parse this AIF.
    if ~(isempty(provided_AIF))
        
        fileID = fopen(provided_AIF, 'r');
        [A, count] = fscanf(fileID, ['%f' ';']);
        fclose(fileID);

        AIF = zeros(tsize,1);
        
        if length(A) > tsize
            err('Provided AIF has more time points than provided DCE signal. Aborting processing');
        elseif length(A) < tsize 
            % Some programs provide AIFs shorter than the input signal.
            % One can either pad in front or in the back. We pad in back by
            % default.

            % Front padding with zeros.
            % AIF((end-length(A)+1):end) = A;
        
            % Back-padding with mean.
            AIF(1:length(A)) = A;
            AIF((length(A)+1):end) = mean(AIF((length(A)-10):length(A)));         
        else
            AIF = A;
        end
        
        if convert_AIF == 0
            gd_AIF = AIF;
        else
            baselineAIF = mean(AIF(first_baseline:last_baseline));
            relAIF = AIF./baselineAIF;
            R1pre = 1./T1_blood;   %1/msec
            a = exp(-TR.*R1pre);
            TERM = (1-a)./(1-a.*cos(alpha_rad)); 
            relAIF = relAIF.*TERM;
            gd_log_term = (relAIF-1)./(a.*(relAIF.*cos(alpha_rad)-1));
            gd_AIF = -(1./(relaxivity*TR)) .* log(gd_log_term);
            gd_AIF = gd_AIF./(1-hematocrit);            
        end

        % Optional simple AIF smoothing. Unsure of usefullness/efficacy at present.
        % gd_AIF = wden(relAIF,'heursure','s','one',5,'sym8');
        
    else
        % If no AIF provided, create a population AIF using the Parker
        % model.
        gd_AIF=generateAIF(tsize,time_interval,bolus_start);
        
        % The Parker model was meant for Gd doses of 1 mmol. We have
        % here an ad-hoc method for other doses, but it isn't perfectly
        % accurate. More precise modifications of the Parker model should
        % be added in the future.
        if (gd_dose ~= 1)
            dose_scalar = 1/gd_dose;
            AIF=generateAIF(int64(tsize*dose_scalar),time_interval,bolus_start);
            newAIF = AIF;
            for p = int64((bolus_start):dose_scalar:(dose_scalar*tsize-(bolus_start)))
                newAIF(bolus_start + (p-bolus_start)/2) = (AIF(p) + AIF(p+1))/2;
            end
            newAIF = newAIF(1:tsize);
            newAIF = newAIF./1.5;
            gd_AIF = newAIF;
        end
    end
 
    % Apply gaussian blurring if specified. Meant to only be applied on
    % axial slices, but so far not very specific.
    if (gaussian_kernel_blur > 0)
        for t = 1:tsize
            for z = 1:zsize
                dce4D.img(:,:,z,t) = imgaussfilt(dce4D.img(:,:,z,t), gaussian_kernel_blur);
            end
        end
    end
    
    % Optionally output the principal components of a PCA analysis.
    if (PCA_output > 0)
    
        PCA_img = dce4D.img;
        PCA_img(isnan(PCA_img)) = 0;
        R = PCA_output;
        data = reshape(double(PCA_img), [xsize*ysize*zsize, tsize]);
        [V, D] = eig(data'*data);
        [D, I] = sort(diag(D), 'descend');
        V = V(:, I);
        U = data*V;
        newdata = reshape (U(:, 1:R), [xsize, ysize, zsize, R]);

        tmp=dce4D;
        tmp.hdr.dime.dim(1)=4;
        tmp.hdr.dime.dim(5)=10;
        tmp.hdr.dime.pixdim(1)=1;
        tmp.hdr.dime.datatype=16;
        tmp.hdr.dime.bitpix=32;  % make sure it is a float image
        tmp.hdr.dime.cal_max=0;
        tmp.hdr.dime.glmax=0;     
        tmp.img=newdata;
        
        fn=['eigs.nii.gz'];
        save_untouch_nii(tmp,strcat(output_path, fn));
        
        % Optionally threshold data values by the second component of the
        % PCA analysis. More testing needs to be done to see if this is a
        % consistent strategy.
        if noise_threshold > 0
            threshdata = newdata(:,:,:,2) > 0;
            testmean = mean2(nonzeros(reshape(newdata(threshdata),1,[]))) * noise_threshold;
        end
    end
    
    % Start parallel pool, and make sure an instance of parallel pool is
    % not already running. If an old instance has the wrong number of
    % workers, restart it.
    if processes > 0
        current_pool = gcp('nocreate');
        if isempty(current_pool)
            parpool(processes)
        elseif current_pool.NumWorkers ~= processes
            delete(gcp('nocreate'))
            parpool(processes)
        end
    end
            
            
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % calculate ktrans, ve at each voxel
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          parfor x=1:xsize
                            for y=1:ysize
                              for z=1:zsize
%                           for x=2
%                             for y=11
%                                 for z=15
%                                 tic
%                                 s=dce4D.img(x,y,:);
%                                 input_signal_4D=zeros(1,1,tsize);
%                                 for n=1:tsize
%                                         input_signal_4D(1,1,n)=s(1,1,n);
%                                 end
%                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 % convert 4D signal to 4D concentration (reference to
%                                 % Step2a_DRO_signal_to_concentration.m)
%                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 signal_4D = input_signal_4D;
%                                 baselineVol_3D = mean(signal_4D(:,:,first_baseline:last_baseline),3);
%                                 if (baselineVol_3D == 0)
%                                 ktransmap(x,y)=0;
%                                 vemap(x,y)=0;
%                                 aucmap(x,y)=0;
%                                 continue
%                                 end
%                                 R1pre = 1./T1_tissue;   %1/msec
%                                 a = exp(-TR.*R1pre);
%                                 TERM = (1-a)./(1-a.*cos(alpha_rad));
%                                 relSignal_4D = zeros(size(signal_4D));
%                                 for i=1:size(signal_4D,3)
%                                     relSignal_4D(:,:,i) = signal_4D(:,:,i)./baselineVol_3D;
%                                 end
%                                 y_4D = relSignal_4D.*(repmat(TERM,[size(signal_4D,1),size(signal_4D,2),size(signal_4D,3)]));
%                                 % Use y_4D to calculate CA concentration:
%                                 gd_conc_4D = zeros(size(relSignal_4D));
%                                 
%                                 for i=1:size(relSignal_4D,3);
%                                     y_3D = squeeze(y_4D(:,:,i));
%                                     gd_log_term = (y_3D-1)./(a.*(y_3D.*cos(alpha_rad)-1));
%                                     gd_conc_4D(:,:,i) = -(1./(relaxivity*TR)) .* log(gd_log_term);
%                                 end
%                                
%                                 
%                                 gd_conc_4D;
%                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 % fit ktrans and Ve (simplex)
%                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 obs_conc=squeeze(gd_conc_4D);
%                                 init_params=[-2, .1];
%                                 fprintf('%d,%d', x,y)
%                                 [kin_par est_conc] = Step4b_Simplex_Fit_Kinetic_Tofts(obs_conc, gd_AIF, init_params, time_interval_mins);
%                                 k1=exp(kin_par(1)); ktrans=k1;
%                                 Ve=1/(1+exp(-kin_par(2)));
%                                 k2 = k1/Ve;
%                                 init_params=[kin_par(1) kin_par(2)];
% %                                 toc
%                                 if (ktrans > 5)
%                                     ktrans = 0;
%                                 end
%                                 if (Ve > 5)
%                                     Ve = 0;
%                                 end
%                                 fprintf('at (%d, %d), Ve=%f, ktrans=%f\n', x, y, Ve, ktrans);
%                                 ktransmap(x,y)=ktrans;
%                                 vemap(x,y)=Ve;
%                                 aucmap(x,y)=trapz(obs_conc)/trapz(AIF);
%                                 if ~(isempty(mask_file))
%                                     if (mask_img(x,y,z) == 0)
%                                         continue
%                                     end
%                                 end

                                s=dce4D.img(x,y,z,:);

                                if newdata(x,y,z,2) >= testmean
%                                     dce4D.img(x,y,z,:) = mean([dce4D.img(abs(x+1),y,z,:), dce4D.img(abs(x-1),y,z,:), dce4D.img(x,abs(y+1),z,:), dce4D.img(x,abs(y-1),z,:)]);
                                        ktransmap(x,y,z)=-.01;
                                        vemap(x,y,z)=-.01;
                                    continue
                                end
                                
                                if any(isnan(s))
                                    continue
                                end
                                
                                input_signal_4D=zeros(1,1,1,tsize);
                                for n=1:tsize
                                        input_signal_4D(1,1,1,n)=s(1,1,1,n);
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % convert 4D signal to 4D concentration (reference to
                                % Step2a_DRO_signal_to_concentration.m)
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                signal_4D = input_signal_4D;
                                baselineVol_3D = mean(signal_4D(:,:,:,first_baseline:last_baseline),4);
                                
%                                 img_slice = dce4D.img(x,:,:,1);
%                                 meanthresh = mean(img_slice(:));
                                

                                if (baselineVol_3D < 10)
                                ktransmap(x,y,z)=0;
                                vemap(x,y,z)=0;
                                aucmap(x,y,z)=0;
                                continue
                                end
                                
                                
                                a = exp(-TR.*R1_pre(x,y,z));
                                TERM = (1-a)./(1-a.*cos(alpha_rad));
                                relSignal_4D = zeros(size(signal_4D));
                                for i=1:size(signal_4D,4)
                                    relSignal_4D(:,:,:,i) = signal_4D(:,:,:,i)./baselineVol_3D;
                                end
                                y_4D = relSignal_4D.*(repmat(TERM,[size(signal_4D,1),size(signal_4D,2),size(signal_4D,3),size(signal_4D,4)]));
                                % Use y_4D to calculate CA concentration:
                                gd_conc_4D = zeros(size(relSignal_4D));
                                
                                for i=1:size(relSignal_4D,4);
                                    y_3D = squeeze(y_4D(:,:,:,i));
                                    gd_log_term = (y_3D-1)./(a.*(y_3D.*cos(alpha_rad)-1));
                                    gd_conc_4D(:,:,:,i) = -(1./(relaxivity*TR)) .* log(gd_log_term);
                                end

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % fit ktrans and Ve (simplex)
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                obs_conc=squeeze(gd_conc_4D);
                                init_params=[1, 1];
                                [kin_par est_conc] = Step4b_Simplex_Fit_Kinetic_Tofts(obs_conc, gd_AIF, init_params, time_interval_mins);
                                k1=exp(kin_par(1)); ktrans=k1;
                                Ve=1/(1+exp(-kin_par(2)));
                                k2 = k1/Ve;
                                init_params=[kin_par(1) kin_par(2)];
%                                 toc
                                if (ktrans > 3)
                                    ktrans = -.01;
                                end
                                if (Ve > 1)
                                    Ve = 1;
                                end
                                fprintf('at (%d, %d, %d), Ve=%f, ktrans=%f\n', x, y,z, Ve, ktrans);
                                ktransmap(x,y,z)=ktrans;
                                vemap(x,y,z)=Ve;
%                                 aucmap(x,y,z)=trapz(obs_conc(lpbs:lpbs+10))/trapz(gd_AIF(lpbs:lpbs+10));

                              end
                            end
                          end
                else
                    fprintf('File path cannot be read.')
end
    
                        % save ktrans, ve, auc maps
                        tmp=dce4D;
                        tmp.hdr.dime.dim(1)=3;
                        tmp.hdr.dime.dim(5)=1;
                        tmp.hdr.dime.pixdim(1)=1;
                        tmp.hdr.dime.datatype=16;
                        tmp.hdr.dime.bitpix=32;  % make sure it is a float image
                        tmp.hdr.dime.cal_max=0;
                        tmp.hdr.dime.glmax=0;
                        toc
                                    
                        tmp.img=ktransmap;
                        fn=['ktrans.nii.gz']
                        save_untouch_nii(tmp,strcat(output_path, fn));


                        tmp.img=vemap;
                        fn=['ve.nii.gz']
                        save_untouch_nii(tmp,strcat(output_path, fn));

%                         tmp.img=aucmap;
%                         fn=['auc.nii.gz']
%                         save_untouch_nii(tmp,strcat(output_path, fn));
                        
%                         tmp.img = tmp.img(:,:,1)
%                         fn=['dce_01.nii.gz']
%                         save_untouch_nii(tmp,strcat(output_path, fn));
                                     
		end                     
