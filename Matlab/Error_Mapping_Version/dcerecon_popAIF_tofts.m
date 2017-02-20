function dcerecon_popAIF(input4Dfile, outputpath, maskfile, T1Map)

if (nargin==2)
	maskfile='';
end

% load dependent parameters
T1_tissue = 1000; %assumed T1 in tissue
T1_fixed = T1_tissue; %in milliseconds - NOTE THAT THIS MAY BE USED FOR THE CONCENTRATION CALCULATION INSTEAD OF T1 MAPS
r1_Cagent = 0.0045; %Relaxometry of Contrast Agent used jkc
T1_blood = 1440;
hematocrit = .45;
TR  = 5; %5; % Repetition time in millisec!
alpha=30;alpha_rad = (pi/180).*alpha;
Total_scan_time_mins = 11; %11; %jkc
scan_start_time_seconds = 60;
firstbaseline = 1;
firstBaseline = firstbaseline;
init_params=[1 1];


                if (exist(input4Dfile))
                        dce4D=load_untouch_nii(input4Dfile);
                        dce4D.img=dce4D.img(:,11:80,:);
                        size(dce4D.img)
                        xsize=size(dce4D.img,1);
                        ysize=size(dce4D.img,2);
                        tsize=size(dce4D.img,3);
                    duration_seconds = 60*Total_scan_time_mins;
                    FR = duration_seconds/(tsize);

                    mapping = 0;
                    if ~(strcmp('', T1Map))
                       T1Map_file = load_untouch_nii(T1Map);
                       T1Map_img = T1Map_file.img(:,:,:);
                       mean(T1Map_img, 1);
                       mean(dce4D.img(:,:,end), 1);
                       size(T1Map_img);
                       size(dce4D.img);
                       mapping = 1;
%                      T1Map_img(:,:) = imgaussfilt(T1Map_img(:,:), 1);

                        
                        R1_pre = 1./T1Map_img;
                    end
                    if (mapping == 0)
                        R1_pre = (1./T1_fixed).*ones([xsize ysize]);
                    end

                    
                    lastbaseline = round(scan_start_time_seconds/FR);
                    lastBaseline = lastbaseline;
                    lpbs=lastbaseline+1;     %dynamic #
                    FR_mins = FR/60;         % in minutes

                        error_interval = 50;
                        error_start = 510;
                        error_end = 2500;
                        error_test_range = error_start:error_interval:error_end;
                        errormap=zeros(xsize,ysize,size(error_test_range,2));
                        ktransmap=zeros(xsize,ysize);
                        vemap=zeros(xsize,ysize);
                        aucmap=zeros(xsize,ysize);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % generate aif (popAIF)
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
                            
                        dce4D.img = double(dce4D.img);
                        
%                         AIF=generateAIF(tsize,FR,lpbs);
%                         gd_AIF = AIF;

                        AIF = dce4D.img(:,61:70,:);
                        [AIF]=mean(AIF(:,:,:),1);
                        [AIF]=mean(AIF(:,:,:),2);

                        lev = 5;
%                         AIF(1,1,:) = wden(AIF(1,1,:),'heursure','s','one',lev,'sym8');
                        [AIF]=double(AIF(:));

%                         sigma = 700;
%                         AIF = AIF + sigma*randn(size(AIF));
%                         plot(AIF)

                        baselineAIF = mean(AIF(firstBaseline:lastBaseline));
                        relAIF = AIF./baselineAIF;
                        
                        real_relAIF = relAIF;
                        
                        if mapping == 1
                        T1_AIF = T1Map_img(:,1:10);
                        T1_blood=mean(T1_AIF);
                        T1_blood=mean(T1_blood);
                        end
                        
                        R1pre = 1./T1_blood;   %1/msec

                        a = exp(-TR.*R1pre);
                        TERM = (1-a)./(1-a.*cos(alpha_rad)); 
                        relAIF = relAIF.*TERM;
                        
                        gd_log_term = (relAIF-1)./(a.*(relAIF.*cos(alpha_rad)-1));
                        gd_AIF = -(1./(r1_Cagent*TR)) .* log(gd_log_term);
                        gd_AIF = gd_AIF./(1-hematocrit);

%                                                 for t = 1:tsize
%                            dce4D.img(:,:,t) = imgaussfilt(dce4D.img(:,:,t), 1);
%                                                 end
                      
                        
%                         PCA Decomposition
%                         R = 3; % Number of components to keep
%                         data = reshape(double(dce4D.img), [xsize*ysize, tsize]);
%                         [V, D] = eig(data'*data);
%                         [D, I] = sort(diag(D), 'descend');
%                         D = diag(D);
%                         V = V(:, I);
%                         U = data*V;
%                         newdata = reshape (U(:, 1:R), [xsize, ysize, R]);
%                         newdata = reshape(newdata, [], R);
%                         data = newdata*V(:,1:R)';
%                         data = reshape(data, [xsize, ysize, tsize]);
%                         dce4D.img = data;
                        
%                         error_interval = 50;
%                         error_start = 500;
%                         error_end = 2500;
                        error_interval = .01;
                        error_start = .3;
                        error_end = 1;                        
                        error_test_range = error_start:error_interval:error_end;
                        errormap=zeros(xsize,ysize,size(error_test_range,2));
                        real_AIF = gd_AIF; 
                        
                        
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         calculate ktrans, ve at each voxel
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
%                           gd_AIF(142:170) = .7*real_AIF(142:170);

%                           parfor x=1:xsize
%                             for y=1:ysize
                          for x=5
                            for y=45
                              for error_param = error_test_range
                                error_param

                                
%                                 R1_pre = (1./error_param).*ones([xsize ysize]);
                                gd_AIF(142:170) = error_param*real_AIF(142:170);
tic
                                s=double(dce4D.img(x,y,:));

                                input_signal_4D=zeros(1,1,tsize);
                                for n=1:tsize
                                        input_signal_4D(1,1,n)=s(1,1,n);
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % convert 4D signal to 4D concentration (reference to
                                % Step2a_DRO_signal_to_concentration.m)
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                signal_4D = input_signal_4D;
                                baselineVol_3D = mean(signal_4D(:,:,firstBaseline:lastBaseline),3);
                                if (baselineVol_3D == 0)
                                ktransmap(x,y)=0;
                                vemap(x,y)=0;
                                aucmap(x,y)=0;
                                continue
                                end
                                a = exp(-TR.*R1_pre(x,y));
                                TERM = (1-a)./(1-a.*cos(alpha_rad));
                                relSignal_4D = zeros(size(signal_4D));
                                for i=1:size(signal_4D,3)
                                    relSignal_4D(:,:,i) = signal_4D(:,:,i)./baselineVol_3D;
                                end
                                y_4D = relSignal_4D.*(repmat(TERM,[size(signal_4D,1),size(signal_4D,2),size(signal_4D,3)]));
                                % Use y_4D to calculate CA concentration:
                                gd_conc_4D = zeros(size(relSignal_4D));
                                
                                for i=1:size(relSignal_4D,3);
                                    y_3D = squeeze(y_4D(:,:,i));
                                    gd_log_term = (y_3D-1)./(a.*(y_3D.*cos(alpha_rad)-1));
                                    gd_conc_4D(:,:,i) = -(1./(r1_Cagent*TR)) .* log(gd_log_term);
                                end
                                
                                gd_conc_4D;
                                obs_conc=squeeze(gd_conc_4D);
%                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 % fit ktrans and Ve (simplex)
%                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 obs_conc;
%                                 init_params=[-2, 1];
%                                 [kin_par est_conc] = Step4b_Simplex_Fit_Kinetic_Tofts(obs_conc, gd_AIF, init_params, FR_mins);
%                                 k1=exp(kin_par(1)); ktrans=k1;
%                                 Ve=1/(1+exp(-kin_par(2)));
%                                 k2 = k1/Ve;
%                                 init_params=[kin_par(1) kin_par(2)];
% %                                 toc
%                                 if (ktrans > .9)
%                                     ktrans = 0;
%                                 end
%                                 if (Ve > 1)
%                                     Ve = 0;
%                                 end
%                                 fprintf('at (%d, %d), Ve=%f, ktrans=%f\n', x, y, Ve, ktrans);
%                                 ktransmap(x,y)=ktrans;
%                                 vemap(x,y)=Ve;
% %                                 aucmap(x,y)=trapz(obs_conc)/trapz(AIF);                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % characterize error map
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                t = 0:(length(obs_conc)-1);
                                t = t.*(FR_mins);
                                t = t(2);
                                ve_range = .01:(.5/50):.5;
                                ktrans_range = .01:(.35/60):.355;
                                AIF = gd_AIF;
                                for ve_param = ve_range
                                    for ktrans_param = ktrans_range
                                            ktrans = ktrans_param;
                                            ve = ve_param;
                                            kep=ktrans/ve;
                                            est_conc=zeros(1,length(AIF));
                                            for i=2:size(AIF,1)
                                               log_e=-kep*t;
                                               e=exp(log_e);
                                               terma=AIF(i)*(e-log_e-1);
                                               termb=AIF(i-1)*(e-(e*log_e)-1);
                                               integral=(terma-termb)/power(log_e,2);
                                               est_conc(i)=est_conc(i-1)*e + ktrans*t*integral;
                                            end

                                            if size(est_conc,2)~=1
                                                est_conc=est_conc';
                                            end
                                            d = (obs_conc-est_conc).^2;
                                            cost = sum(d);
%                                             if (ve_param == .05 && ktrans_param > .2)
%                                             size(obs_conc)
%                                             size(est_conc)
%                                             timeseries = linspace(0,11,length(obs_conc));
%                                             plot(timeseries, obs_conc, timeseries, est_conc);
%                                             ve_param
%                                             ktrans_param
%                                             a = 1;
%                                             fd = gf
%                                             end
                                            if (cost>1)
                                                cost = 1;
                                            end
                                        size(ktransmap);
                                    	ktransmap(int8((ve-.01)/(.5/50)+1),int8((ktrans-.01)/(.35/60)+1))=cost;
                                        vemap(int8((ve-.01)/(.5/50)+1),int8((ktrans-.01)/(.35/60)+1))=cost;
                                   end                                    
                                end
                                errormap(:,:,int8((error_param-error_start+error_interval)/error_interval)) = ktransmap;
                            end
                          end

                          end

                        assignin('base','ktransmap',ktransmap(:,1:60))
                        assignin('base','vemap',vemap(:,1:60))
%                         contourf(.01:(.5/50):.5,fliplr(.01:(.35/60):.355),flipud(transpose(ktransmap(:,1:60))), 100)              
                        toc

                        
                        % save ktrans, ve, auc maps
                        tmp=dce4D;
                        tmp.hdr.dime.dim(1)=3;
                        tmp.hdr.dime.dim(2)=size(errormap,1);
                        tmp.hdr.dime.dim(3)=size(errormap,2);                     
                        tmp.hdr.dime.dim(4)=size(errormap,3);
                        tmp.hrd.dime.dim(5)=1;
%                         tmp.hdr.dime.dim(1)=2;
%                         tmp.hdr.dime.dim(2)=size(ktransmap,1);
%                         tmp.hdr.dime.dim(3)=size(ktransmap,2);                     
%                         tmp.hdr.dime.dim(4)=1;
%                         tmp.hrd.dime.dim(5)=1;                        
                        tmp.hdr.dime.pixdim(1)=1;
                        tmp.hdr.dime.datatype=16;
                        tmp.hdr.dime.bitpix=32;  % make sure it is a float image
                        tmp.hdr.dime.cal_max=0;
                        tmp.hdr.dime.glmax=0;
            
                        path = input4Dfile;

                        tmp.img=errormap;
                        fn=['errormap.nii.gz']
                        fn = [strcat(outputpath, fn)]
                        save_untouch_nii(tmp,fn);
                        
%                         tmp.img=ktransmap;
%                         fn=['ktrans.nii.gz']
%                         fn = [strcat(outputpath, fn)]
%                         save_untouch_nii(tmp,fn);
% %                         fn=['ktrans.mat']
% %                         fn = [strcat(outputpath, fn)]
% % %                         dicomwrite(tmp.img, fn, 'CreateMode','Copy')
% %                         save(fn, 'ktransmap')
% 
%                         tmp.img=vemap;
%                         fn=['ve.nii.gz']
%                         save_untouch_nii(tmp,strcat(outputpath, fn));
%                         fn=['ve.mat']
%                         fn = [strcat(outputpath, fn)]
% %                         dicomwrite(tmp.img, fn)
%                          
%                         save(fn, 'vemap')   
%                         tmp.img=aucmap;
%                         fn=['auc.nii.gz']
%                         save_untouch_nii(tmp,strcat(outputpath, fn));
                                     
		end                     
