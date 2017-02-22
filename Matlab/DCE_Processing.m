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
    time_interval_mins = time_interval_seconds/60;
    last_baseline = round(bolus_arrival_time_seconds/time_interval_seconds);
    bolus_start = last_baseline + 1;

    % note: make this work irrespective of dimensions
    ktransmap=zeros(xsize,ysize,zsize);
    vemap=zeros(xsize,ysize,zsize);
    aucmap=zeros(xsize,ysize,zsize);

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
        gd_AIF=generateAIF(tsize,time_interval_seconds,bolus_start);
        
        % The Parker model was meant for Gd doses of 1 mmol. We have
        % here an ad-hoc method for other doses, but it isn't perfectly
        % accurate. More precise modifications of the Parker model should
        % be added in the future.
        if (gd_dose ~= 1)
            dose_scalar = 1/gd_dose;
            AIF=generateAIF(int64(tsize*dose_scalar),time_interval_seconds,bolus_start);
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
    
    % Convert input signal volume
    baselineVol_3D = mean(dce4D.img(:,:,:,first_baseline:last_baseline),4);
    a = exp(-TR.*R1_pre);
    TERM = (1-a)./(1-a.*cos(alpha_rad));
    relSignal_4D = zeros(size(dce4D.img));
    for i=1:size(dce4D.img, 4)
        relSignal_4D(:,:,:,i) = double(dce4D.img(:,:,:,i)) ./ baselineVol_3D .* TERM;
        relSignal_4D(:,:,:,i) = (relSignal_4D(:,:,:,i) - 1) ./ (a .* (relSignal_4D(:,:,:,i) .* cos(alpha_rad) - 1));
    end
    relSignal_4D = -(1./(relaxivity*TR)) .* log(relSignal_4D);
    
    % Start parallel pool, and make sure an instance of parallel pool is
    % not already running. If an old instance has the wrong number of
    % workers, restart it.
    if processes > 0
        current_pool = gcp('nocreate');
        if ~(isempty(current_pool))
            if ~(current_pool.NumWorkers == processes)
                delete(gcp('nocreate'));
                parpool(processes);
            end
        else
            parpool(processes);
        end
    end
    
    % note: would be nice for this loop to run regardless of dimension.
    parfor x=1:xsize
        for y=1:ysize
            for z=1:zsize
                
                % Extract time signal
                signal_4D=double(relSignal_4D(x,y,z,:));
                
                % Optional PCA Thresholding
%                 if newdata(x,y,z,2) >= 10
%                 %       dce4D.img(x,y,z,:) = mean([dce4D.img(abs(x+1),y,z,:), dce4D.img(abs(x-1),y,z,:), dce4D.img(x,abs(y+1),z,:), dce4D.img(x,abs(y-1),z,:)]);
%                         ktransmap(x,y,z)=-.01;
%                         vemap(x,y,z)=-.01;
%                     continue
%                 end

                % Ignore NaN Values
                if any(isnan(signal_4D))
                    continue
                end

                % Low intensity threshold.
                % Note: this threshold is arbitrary. It may be better to
                % make it a parameter or automatically set it with regards
                % to the distribution of values within the volume.
                if baselineVol_3D(x,y,z) < 10
                    continue
                end                
                
                % Call on the fitting function to generate ktrans and Ve
                % values.
                observed_concentration=squeeze(signal_4D);
                initial_params=[1, 1];
                [kinetic_params, estimated_concentration] = Simplex_Fit_Kinetic_Tofts(observed_concentration, gd_AIF, initial_params, time_interval_mins);
                ktrans=exp(kinetic_params(1));
                Ve=1/(1+exp(-kinetic_params(2)));
                
                % Warm-Fitting for subsequent parameters. Optional, may
                % increase speed but decrease accuracy. TODO: Improve
                % warm-fitting.
                % initial_params=[kinetic_params(1) kinetic_params(2)];
                
                % Mask degenerate results. These threshold are arbitrary.
                % In particular for ktrans, it may be better 
                if (ktrans > 3)
                    ktrans = -.01;
                end
                if (Ve > 1)
                    Ve = 1;
                end
                
                % TODO: Add options for non-verbose mode. Code runs about
                % 25% faster without printing output.
                fprintf('at (%d, %d, %d), Ve=%f, ktrans=%f\n', x, y,z, Ve, ktrans);
                
                % Save out results.
                ktransmap(x,y,z)=ktrans;
                vemap(x,y,z)=Ve;
                aucmap(x,y,z)=trapz(observed_concentration(last_baseline:tsize))/trapz(gd_AIF(last_baseline:tsize));
            end
        end
    end

end
    
% Save images.
% note: make outputs optionally specified
% and make dimension options below more flexible.
tmp.hdr.dime.dim(1)=3;
tmp.hdr.dime.dim(5)=1;
tmp.hdr.dime.pixdim(1)=1;
tmp.hdr.dime.datatype=16;
tmp.hdr.dime.bitpix=64;  % make sure it is a float image
tmp.hdr.dime.cal_max=0;
tmp.hdr.dime.glmax=0;

tmp.img=ktransmap;
fn=['ktrans.nii.gz']
save_untouch_nii(tmp,strcat(output_path, fn));

tmp.img=vemap;
fn=['ve.nii.gz']
save_untouch_nii(tmp,strcat(output_path, fn));

tmp.img=aucmap;
fn=['auc.nii.gz']
save_untouch_nii(tmp,strcat(output_path, fn));
                                     
end                     
