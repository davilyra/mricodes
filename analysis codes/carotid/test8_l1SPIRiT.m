%% Adding required toolboxes to the path
addpath(genpath('/home/mri-vh/Davi/SPIRiT_v0.1'))
% addpath('/home/mri-vh/Davi/SPIRiT_v0.1/utils')

clear all;
close all;
clc;

% datasets = 12;
datasets = 10:18;
dts = length(datasets);

% data and recon information
ncoils = 15;        % number of coils
N = [128, 256];     % imagesize
nslices = 18;       % number of slices
accelfactor = 8;
% mkdir(['accelfactor_',num2str(accelfactor)]);


%% Sampling
% pseudorandom sampling
load(['samplingmasks/sampling_mask_',num2str(accelfactor),'.mat']);


%% SPIRiT parameters
kSize = [5,5];      % SPIRiT kernel size
nIter = 05;         % number of iteration; phantom requires twice as much as the brain
CalibTyk = 0.1;     % Tykhonov regularization in the calibration
wavWeight = 0.15;   % Wavelet soft-thresholding regularization in the reconstruction


%% SPIRiT Calibration
[CalibSize, dcomp] = getCalibSize(mask);  % get size of calibration area from mask
compmatrix = repmat(dcomp,[1,1,ncoils]);  % compensation matrix


%% Initializing reconstruction matrices
rdc1 = zeros(N(1), N(2), ncoils); % SPIRiT and CS data matrix
rdc3 = zeros(N(1),N(2),nslices);


%% Recon
for kk = 1:dts
    disp(['reconstructing dataset number ',num2str(datasets(kk))]),
    load(['Run2730.6904.',num2str(datasets(kk)),'_2.mat']);
    
    for slice = 1:nslices;
        disp(['working with slice number ',num2str(slice)])
        for i = 1:ncoils
            % Organizing the data
            rd1a = raws2(slice,:,:,i);
            rd1a = squeeze(rd1a);
            rd2a = flipud(rd1a(1:floor((size(rd1a,1) + 1)/2),:));
            rd3a = (rd1a(floor((size(rd1a,1) + 1)/2 + 1):size(rd1a,1),:));
            rda = [rd2a; rd3a];
            rda = rda.*mask;
            
            rd1b = raws2(slice + 18,:,:,i);
            rd1b = squeeze(rd1b);
            rd2b = flipud(rd1b(1:floor((size(rd1b,1) + 1)/2),:));
            rd3b = (rd1b(floor((size(rd1b,1) + 1)/2 + 1):size(rd1b,1),:));
            rdb = [rd2b; rd3b];
            rdb = rdb.*mask;
            
            rdiff = (rda - rdb)./pdf;
            rdc1(:,:,i) = rdiff; % SPIRiT data
        end
        
        DATAcomp = rdc1.*compmatrix;
        scale_fctr = norm(DATAcomp(:))/sqrt(ncoils)/20;
        rdc1 = rdc1./scale_fctr; % SPIRiT and CS Recon
        
        % SPIRiT Recon
        disp('performing calibration for SPIRiT')
        kCalib = crop(rdc1,[CalibSize, ncoils]);
        kernel = zeros([kSize, ncoils, ncoils]);
        
        [AtA,] = corrMatrix(kCalib,kSize);
        for n=1:ncoils
            kernel(:,:,:,n) = calibrate(AtA,kSize,ncoils,n,CalibTyk);
        end
        GOP = SPIRiT(kernel, 'conv',N);
        
        disp('performing SPIRiT reconstruction')
        [res_pocs] = pocsSPIRiT(rdc1,GOP,nIter,rdc1,wavWeight,0);
        im_pocsspirit = ifft2c(res_pocs);
        rdc3(:,:,slice) = sos(im_pocsspirit);
        disp(' ')
        
    end
    
    mip = max(rdc3, [], 3);
    figure; imagesc(abs(mip)),colormap(gray)
    pause(.1)
    
    disp('Done!!')
    save(['accelfactor_',num2str(accelfactor),'/carotid_dataset_',num2str(datasets(kk)),'_l1spirit_recon.mat'], 'rdc3', 'N','nslices');
    disp(' ')

end
