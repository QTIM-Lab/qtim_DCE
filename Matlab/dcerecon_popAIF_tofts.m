function dcerecon_popAIF(input4Dfile, outputpath, maskfile, provided_AIF, T1Map, kernel_blur, PCA, dose, Flip_Angle, Static_T1, Threshold)

tic
if (nargin==2)
	maskfile='';
    provided_AIF='';
elseif (nargin==3)
    maskfile='';
end

suppress_PCA = 10;

% Change these if you are using an ROI for the AIF.
T1_tissue = Static_T1; %assumed T1 in tissue
T1_fixed = T1_tissue; %in milliseconds - NOTE THAT THIS MAY BE USED FOR THE CONCENTRATION CALCULATION INSTEAD OF T1 MAPS
r1_Cagent = 0.0039; %.0037 for competition %.0039 for MGH %Relaxometry of Contrast Agent used jkc
T1_blood = 1440;
hematocrit = .45;

% Change these every time.
TR  = 5.504; %Competition 3.594 %TMZ 5.504 % CED/NHX 6.8 5; % Repetition time in millisec!
alpha=Flip_Angle;alpha_rad = (pi/180).*alpha; %Competition 15 % CED/NHX 10
Total_scan_time_mins = 6.86; %Competition 5 %TMZ 6.866 % CED/NHX 6 %11;%jkc
scan_start_time_seconds = 60; %Competition 40 %TMZ 60 %CED/NHX 160

input4Dfile

firstbaseline = 1;
firstBaseline = firstbaseline;
init_params=[1 1];

R1_pre = 0;
mapping = 0;

if ~(isempty(maskfile))
    mask = load_untouch_nii(maskfile);
    mask_img = mask.img;
end

                if (exist(input4Dfile))
                        dce4D=load_untouch_nii(input4Dfile);
                        size(dce4D.img)
                        xsize=size(dce4D.img,1);
                        ysize=size(dce4D.img,2);
                        zsize=size(dce4D.img,3);
                        tsize=size(dce4D.img,4);
                        
                   if ~(strcmp('', T1Map))
                       T1Map_file = load_untouch_nii(T1Map);
                       T1Map_img = T1Map_file.img;
                       T1Map_img = flip(T1Map_img,2);
                       mapping = 1;

                        if (kernel_blur > 0)
                                            for z = 1:size(T1Map_img,3)
                        T1Map_img(:,:,z) = imgaussfilt(T1Map_img(:,:,z), kernel_blur);
                                            end
                        end

                        R1_pre = 1./T1Map_img;
                   end
                    duration_seconds = 60*Total_scan_time_mins;
                    FR = duration_seconds/tsize;
                    lastbaseline = round(scan_start_time_seconds/FR);
                    lastBaseline = lastbaseline;
                    lpbs=lastbaseline+1;     %dynamic #
                    FR_mins = FR/60;         % in minutes
                    if (mapping == 0)
                        R1_pre = (1./T1_fixed).*ones([xsize ysize zsize]);
                    end
                        ktransmap=zeros(xsize,ysize,tsize);
                        vemap=zeros(xsize,ysize,tsize);
                        aucmap=zeros(xsize,ysize,tsize);
                        
                        ktransmap=zeros(xsize,ysize,zsize,tsize);
                        vemap=zeros(xsize,ysize,zsize,tsize);
                        aucmap=zeros(xsize,ysize,zsize,tsize);
                    
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % generate aif (popAIF)
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if ~(isempty(provided_AIF))
                           fileID = fopen(provided_AIF, 'r');
                           [A, count] = fscanf(fileID, ['%f' ';']);
                           fclose(fileID);
                           
                        A;
                        length(A);
                        AIF = zeros(tsize,1);
%                         AIF((end-length(A)+1):end) = A;
                        AIF(1:length(A)) = A;
                        if length(A) < tsize
                           AIF((length(A)+1):end) = mean(AIF((length(A)-10):length(A)));
                        end
                        
%                         relAIF = AIF;
                        gd_AIF = AIF;
%                         baselineAIF = mean(AIF(firstBaseline:lastBaseline))
%                         relAIF = AIF./baselineAIF;
%                         
%                         R1pre = 1./T1_blood;   %1/msec
%                         a = exp(-TR.*R1pre);
%                         TERM = (1-a)./(1-a.*cos(alpha_rad)); 
%                         relAIF = relAIF.*TERM;
%                         gd_log_term = (relAIF-1)./(a.*(relAIF.*cos(alpha_rad)-1));
%                         gd_AIF = -(1./(r1_Cagent*TR)) .* log(gd_log_term);
%                         gd_AIF = gd_AIF./(1-hematocrit);
                        
                        % gd_AIF = wden(relAIF,'heursure','s','one',5,'sym8');
                        else
                        AIF=generateAIF(tsize,FR,lpbs);
                        gd_AIF = AIF;
                        if (dose<1)
                        AIF=generateAIF(tsize*2,FR,lpbs);
                        newAIF = AIF;
                        for p = lpbs:2:(2*tsize-lpbs)
                            newAIF(lpbs + (p-lpbs)/2) = (AIF(p) + AIF(p+1))/2;
                        end
                        newAIF = newAIF(1:tsize);
                        newAIF = newAIF./1.5;
                        gd_AIF = newAIF;
                        end
                        end
%                         [AIF]=dce4D.img(19,75,:);
%                         [AIF]=double(AIF(:));
%                         
%                         baselineAIF = mean(AIF(firstBaseline:lastBaseline));
%                         relAIF = AIF./baselineAIF;
%                         
%                         R1pre = 1./T1_blood;   %1/msec
%                         a = exp(-TR.*R1pre);
%                         TERM = (1-a)./(1-a.*cos(alpha_rad)); 
%                         relAIF = relAIF.*TERM;
%                         
%                         gd_log_term = (relAIF-1)./(a.*(relAIF.*cos(alpha_rad)-1));
%                         gd_AIF = -(1./(r1_Cagent*TR)) .* log(gd_log_term);
%                         gd_AIF = gd_AIF./(1-hematocrit);
% 
                            if (kernel_blur > 0)
                                                for t = 1:tsize
                                                    for z = 1:zsize
                           dce4D.img(:,:,z,t) = imgaussfilt(dce4D.img(:,:,z,t), kernel_blur);
                                                    end
                                                end
                            end
                         
                          if (PCA > 0)
                        dce4D.img(isnan(dce4D.img)) = 0;    
                        R = PCA; % Number of components to keep
                        data = reshape(double(dce4D.img), [xsize*ysize*zsize, tsize]);
                        [V, D] = eig(data'*data);
                        [D, I] = sort(diag(D), 'descend');
                        D = diag(D);
                        V = V(:, I);
                        U = data*V;
                        newdata = reshape (U(:, 1:R), [xsize, ysize, zsize, R]);
                        newdata = reshape(newdata, [], R);
                        data = newdata*V(:,1:R)';
                        data = reshape(data, [xsize, ysize, zsize, tsize]);
                        dce4D.img = data;                             
                        
                        
                          end
                          
                          if (suppress_PCA > 0)
                             
                              
                        dce4D.img(isnan(dce4D.img)) = 0;
                          R = suppress_PCA;
                        data = reshape(double(dce4D.img), [xsize*ysize*zsize, tsize]);
                        [V, D] = eig(data'*data);
                        [D, I] = sort(diag(D), 'descend');
%                         D
                        V = V(:, I);
                        U = data*V;
                        size(U)
                        newdata = reshape (U(:, 1:R), [xsize, ysize, zsize, R]);
                        size(newdata)
                        
%                         for i = 1:10
%                         mean2(nonzeros(abs(newdata(:,:,:,i))))
%                         size(V)
%                         plot(V(:,i))
%                         pause
%                         end
                                                
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
                        save_untouch_nii(tmp,strcat(outputpath, fn));
%                         newdata = reshape(newdata, [], R);
%                         size(newdata)
%                         data = newdata*V(:,1:R)';
%                         size(data)
%                         data = reshape(data, [xsize, ysize, zsize, tsize]);   
%                         size(data)
                        threshdata = newdata(:,:,:,2) > 0;
                        testmean = mean2(nonzeros(reshape(newdata(threshdata),1,[]))) * Threshold;
%                         fd=gf
                        %timeseries = linspace(0,60,60);
                        %plot(timeseries, squeeze(dce4D.img(70,83,5,:)),timeseries,squeeze(data(70,83,5,:)))
                        %plot()
                        %image(dce4D.img(:,:,5,50))
                        %image(data(:,:,10,40))
                        %fd=gd
                        
                        
                          end
                          
                          fn=['original.nii.gz'];
                          save_untouch_nii(dce4D,strcat(outputpath, fn));
                          
                          gd_AIF;
                          
                          
%                         plot(gd_AIF)
%                           FD = GD
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
%                                 baselineVol_3D = mean(signal_4D(:,:,firstBaseline:lastBaseline),3);
%                                 if (baselineVol_3D == 0)
%                                 ktransmap(x,y)=0;
%                                 vemap(x,y)=0;
%                                 aucmap(x,y)=0;
%                                 continue
%                                 end
%                                 R1pre = 1./T1_fixed;   %1/msec
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
%                                     gd_conc_4D(:,:,i) = -(1./(r1_Cagent*TR)) .* log(gd_log_term);
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
%                                 [kin_par est_conc] = Step4b_Simplex_Fit_Kinetic_Tofts(obs_conc, gd_AIF, init_params, FR_mins);
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
%                                 if ~(isempty(maskfile))
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
                                baselineVol_3D = mean(signal_4D(:,:,:,firstBaseline:lastBaseline),4);
                                
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
                                    gd_conc_4D(:,:,:,i) = -(1./(r1_Cagent*TR)) .* log(gd_log_term);
                                end

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % fit ktrans and Ve (simplex)
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                obs_conc=squeeze(gd_conc_4D);
                                init_params=[1, 1];
                                [kin_par est_conc] = Step4b_Simplex_Fit_Kinetic_Tofts(obs_conc, gd_AIF, init_params, FR_mins);
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
                        save_untouch_nii(tmp,strcat(outputpath, fn));


                        tmp.img=vemap;
                        fn=['ve.nii.gz']
                        save_untouch_nii(tmp,strcat(outputpath, fn));

%                         tmp.img=aucmap;
%                         fn=['auc.nii.gz']
%                         save_untouch_nii(tmp,strcat(outputpath, fn));
                        
%                         tmp.img = tmp.img(:,:,1)
%                         fn=['dce_01.nii.gz']
%                         save_untouch_nii(tmp,strcat(outputpath, fn));
                                     
		end                     
