clc%% Adding required toolboxes to the path
addpath(genpath('/home/mri-vh/Davi/SPIRiT_v0.1'))
addpath(genpath('/home/mri-vh/Davi/functions'))

clear all;
% close all;
clc;

% Data parameters
load('params.mat');
rdc30 = zeros(params(2),params(1),params(3));
% rdc31 = zeros(params(2),params(1),params(3));
% rdc32 = zeros(params(2),params(1),params(3));


% SPIRiT parameters
kSize = [3,3];      % SPIRiT kernel size
nIter = 3;         % number of iteration; phantom requires twice as much as the brain
CalibTyk = 0.1;     % Tykhonov regularization in the calibration
wavWeight = 0.03125;   % Wavelet soft-thresholding regularization in the reconstruction
% wavWeight  = 0;

% recon_data = zeros(params(2),params(1),params(4),params(3),params(5));
recon_data1 = zeros(params(2),params(1),params(4),params(3),params(5));

type = 1; %[1] = pese; [2] = pero

% data and recon information
% accelfactors = [2 4 5 8 10];
accelfactors = [10];
uslth = length(accelfactors);

for iii = 1:uslth
    accelfactor = accelfactors(iii);

switch type
    case 1
        application = 'FBI/pese';
        mkdir([application(5:end),'/accelfactor_',num2str(accelfactor)]);
        
        %% Sampling
        % pseudorandom sampling
        load(['/home/mri-vh/Davi/samplingmasks/',application,'/sampling_mask_',num2str(accelfactor),'.mat']);
        maska = permute(repmat(mask,[1,1,params(2),params(4)]),[2 1 3 4]);
        
        
        %% Loading data
        load([application(5:end),'/accelfactor_',num2str(accelfactor),...
            '/rawdata_fbi1_',application(5:end),'_us_',num2str(accelfactor),'.mat']);
        rd_us_comp = rd_fbi_us; % applying k-space density compensation
%         rdt_mean = mean(rd_us_comp(:,:,:,:,2:25),5);
        rds_mean = mean(rd_us_comp(:,:,:,:,:),3);
        
%         rdt_mean = sum(rd_us_comp(:,:,:,:,2:25),5)/24;
%         rds_mean = sum(rd_us_comp(:,:,:,:,:),3)/10;
        clear rd_fbi_us
        
        
        %% Recon
        for kk = 2:params(5)
%         for kk = 15
            disp(['reconstructing bbti number ',num2str(kk)])
            
            for slice = 1:params(3);
                
                %% SPIRiT Calibration
                [~, dcomp] = getCalibSize(squeeze(permute(maska(slice,:,:,1),[3 2 1 4 5])));  % get size of calibration area from mask
                CalibSize = [5,5];
                compmatrix = repmat(dcomp,[1,1,params(4)]);  % compensation matrix
                
                disp(['working with slice number ',num2str(slice)])
                
                rdc1 = squeeze(rd_us_comp(:,:,slice,:,kk));
                DATAcomp = rdc1.*compmatrix;
                scale_fctr = norm(DATAcomp(:))/sqrt(params(4))/20;
                rdc1 = rdc1./scale_fctr; % SPIRiT and CS Recon
                
%                 rdc0 = squeeze(rdt_mean(:,:,slice,:));
%                 rdc0 = rdc0./scale_fctr;
%                 
                rdc2 = squeeze(rds_mean(:,:,1,:,kk));
                rdc2 = rdc2./scale_fctr;
                
                % SPIRiT Recon
                disp('performing calibration for SPIRiT')
%                 kCalib = crop(rdc1,[CalibSize, params(4)]);
                kCalib = crop(rdc2,[CalibSize, params(4)]);
                kernel = zeros([kSize, params(4), params(4)]);
                
                [AtA,] = corrMatrix(kCalib,kSize);
                for n=1:params(4)
                    kernel(:,:,:,n) = calibrate(AtA,kSize,params(4),n,CalibTyk);
                end
                GOP = SPIRiT(kernel, 'conv', [params(1), params(2)]);
                
                masku = permute(squeeze(maska(slice,:,:,:)),[2 1 3]);
                
                disp('performing SPIRiT reconstruction')
%                 [res_pocs] = pocsSPIRiT(rdc1,GOP,nIter,rdc1,wavWeight,0);
                [res_pocs0] = pocsSPIRiT2(rdc1,GOP,nIter,rdc1,wavWeight,masku,0);
                im_pocsspirit0 = ifft2c(res_pocs0);
                recon_data1(:,:,:,slice,kk) = res_pocs0;
                rdc30(:,:,slice) = combinechannels((im_pocsspirit0),3,2);
                
%                 [res_pocs1] = pocsSPIRiT2(rdc1,GOP,nIter,rdc0,wavWeight,masku,0);
%                 im_pocsspirit1 = ifft2c(res_pocs1);
%                 recon_data(:,:,:,slice,kk) = res_pocs1;
%                 rdc31(:,:,slice) = combinechannels((im_pocsspirit1),3,2);

%                 [res_pocs2] = pocsSPIRiT2(rdc1,GOP,nIter,rdc2,wavWeight,masku,0);
%                 im_pocsspirit2 = ifft2c(res_pocs2);
%                 recon_data(:,:,:,slice,kk) = res_pocs2;
%                 rdc32(:,:,slice) = combinechannels((im_pocsspirit2),3,2);

                disp(' ')
                
            end
            
            mip0 = max(rdc30, [], 3);
%             mip1 = max(rdc31, [], 3);
%             mip2 = max(rdc32, [], 3);
% %             figure;
%             subplot(1,3,1); 
            imshow(abs(mip0),[]), colormap(gray)
% %             subplot(1,3,2);
%             imshow(abs(mip1),[]), colormap(gray)
% %             subplot(1,3,3);
%             imshow(abs(mip2),[]), colormap(gray)

%             subplot(1,3,1); imagesc(abs(mip0)), colormap(gray)
%             subplot(1,3,2); imagesc(abs(mip1)), colormap(gray)
%             subplot(1,3,3); imagesc(abs(mip2)), colormap(gray)
            
            pause(.1)
            
            disp('Done!!')
            disp(' ')
            
        end
%         save([application(5:end),'/accelfactor_',num2str(accelfactor),'/fbi_l1_spirit_recon.mat'], 'recon_data1','-v7.3');
%         save([application(5:end),'/accelfactor_',num2str(accelfactor),'/fbi_l1_wavelets_spirit_recon.mat'], 'recon_data1','-v7.3');

        %         save([application(5:end),'/accelfactor_',num2str(accelfactor),'/fbi_l1_t-spirit_recon.mat'], 'recon_data','-v7.3');        
        save([application(5:end),'/accelfactor_',num2str(accelfactor),'/fbi_s-spirit_recon.mat'], 'recon_data1','-v7.3');
        
        
    case 2
        application = 'FBI/pero';
        mkdir([application(5:end),'/accelfactor_',num2str(accelfactor)]);
        
        
        %% Sampling
        % pseudorandom sampling
        load(['/home/mri-vh/Davi/samplingmasks/',application,'/sampling_mask_',num2str(accelfactor),'.mat']);
        maska = permute(repmat(mask,[1,1,params(3),params(4)]),[2 1 3 4]);
        
        
        %% SPIRiT Calibration
        [CalibSize, dcomp] = getCalibSize(mask);  % get size of calibration area from mask
        compmatrix = permute(repmat(dcomp,[1,1,params(4)]),[2 1 3]);  % compensation matrix
        
        
        %% Loading data
        load([application(5:end),'/accelfactor_',num2str(accelfactor),...
            '/rawdata_fbi1_',application(5:end),'_us_',num2str(accelfactor),'.mat']);
        rd_us_comp = rd_fbi_us; % applying k-space density compensation

%         rdt_mean = mean(rd_us_comp(:,:,:,:,2:25),5);
%         rds_mean = mean(rd_us_comp(:,:,:,:,:),3);
        
%         rdt_mean = sum(rd_us_comp(:,:,:,:,2:25),5)/24;
%         rds_mean = sum(rd_us_comp(:,:,:,:,:),3)/10;
        clear rd_fbi_us
        
        
        %% Recon
        for kk = 2:params(5)
            disp(['reconstructing bbti number ',num2str(kk)])
            
            for slice = 1:params(3);
                disp(['working with slice number ',num2str(slice)])
                
                rdc1 = squeeze(rd_us_comp(:,:,slice,:,kk));
                DATAcomp = rdc1.*compmatrix;
                scale_fctr = norm(DATAcomp(:))/sqrt(params(4))/20;
                rdc1 = rdc1./scale_fctr;
                
%                 rdc0 = squeeze(rdt_mean(:,:,slice,:));
%                 rdc0 = rdc0./scale_fctr;
%                 
%                 rdc2 = squeeze(rds_mean(:,:,1,:,kk));
%                 rdc2 = rdc2./scale_fctr;
                
                
                % SPIRiT Recon
                disp('performing calibration for SPIRiT')
                kCalib = crop(rdc1,[CalibSize, params(4)]);
                kernel = zeros([kSize, params(4), params(4)]);
                
                [AtA,] = corrMatrix(kCalib,kSize);
                for n=1:params(4)
                    kernel(:,:,:,n) = calibrate(AtA,kSize,params(4),n,CalibTyk);
                end
                GOP = SPIRiT(kernel, 'conv', [params(1), params(3)]);
                
                masku = squeeze(maska(:,:,slice,:));
                
                disp('performing SPIRiT reconstruction')
                [res_pocs0] = pocsSPIRiT2(rdc1,GOP,nIter,rdc1,wavWeight,masku,0);
                im_pocsspirit0 = ifft2c(res_pocs0);
                recon_data1(:,:,:,slice,kk) = res_pocs0;
                rdc30(:,:,slice) = combinechannels((im_pocsspirit0),3,2);
                
%                 [res_pocs1] = pocsSPIRiT2(rdc1,GOP,nIter,rdc0,wavWeight,masku,0);                
%                 im_pocsspirit1 = ifft2c(res_pocs1);
% %                 recon_data(:,:,:,slice,kk) = res_pocs1;
%                 rdc31(:,:,slice) = combinechannels((im_pocsspirit1),3,2);
%                 
%                 [res_pocs2] = pocsSPIRiT2(rdc1,GOP,nIter,rdc2,wavWeight,masku,0);
%                 im_pocsspirit2 = ifft2c(res_pocs2);
% %                 recon_data(:,:,:,slice,kk) = res_pocs2;
%                 rdc32(:,:,slice) = combinechannels((im_pocsspirit2),3,2);
                
                
                disp(' ')
                
            end
            
            mip0 = max(rdc30, [], 3);
%             mip1 = max(rdc31, [], 3);
%             mip2 = max(rdc32, [], 3);
%             figure;
%             subplot(1,3,1);
            imshow(abs(mip0),[]), colormap(gray)
            
%             subplot(1,3,2);
%             imshow(abs(mip1),[]), colormap(gray)
%             subplot(1,3,3);
%             imshow(abs(mip2),[]), colormap(gray)
%             subplot(1,3,1); imagesc(abs(mip0)), colormap(gray)
%             subplot(1,3,2); imagesc(abs(mip1)), colormap(gray)
%             subplot(1,3,3); imagesc(abs(mip2)), colormap(gray)
            pause(.1)
            
            disp('Done!!')
            disp(' ')
            
        end
%         save([application(5:end),'/accelfactor_',num2str(accelfactor),'/fbi_l1_spirit_recon.mat'], 'recon_data1','-v7.3');
        save([application(5:end),'/accelfactor_',num2str(accelfactor),'/fbi_l1_wavelets_spirit_recon.mat'], 'recon_data1','-v7.3');
%         save([application(5:end),'/accelfactor_',num2str(accelfactor),'/fbi_t-spirit_recon.mat'], 'recon_data','-v7.3');
%         save([application(5:end),'/accelfactor_',num2str(accelfactor),'/fbi_s-spirit_recon.mat'], 'recon_data','-v7.3');
        
end

end