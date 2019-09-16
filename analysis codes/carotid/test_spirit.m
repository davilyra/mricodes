addpath('/home/mri-vh/Davi/SPIRiT_v0.1')
addpath('/home/mri-vh/Davi/SPIRiT_v0.1/utils')

clear all;
% close all;
clc;

% datasets = 13;
datasets = 10:18;
dts = length(datasets);

% data and recon information
ncoils = 15;        % number of coils
N = [128, 256];     % imagesize
kSize = [7,7];      % SPIRiT kernel size
nIter = 05;         % number of iteration; phantom requires twice as much as the brain
CalibTyk = 0.01;    % Tykhonov regularization in the calibration
wavWeight = 0.015;  % Wavelet soft-thresholding regularization in the reconstruction


% pseudorandom sampling
pdf = genPDF(N, 15, 5/8, 2, 0, 0);
mask = genSampling(pdf,500,2);

[CalibSize, dcomp] = getCalibSize(mask);  % get size of calibration area from mask
compmatrix = repmat(dcomp,[1,1,ncoils]);   % compensation matrix

rdc1 = zeros(N(1), N(2), ncoils); % data matrix
% rdc2 = zeros(N(1), N(2), ncoils); % data matrix

figure; clf;
figr = floor(sqrt(dts));
figc = ceil(dts/figr);

    
for kk = 1:dts
    disp(['reconstructing dataset number ',num2str(datasets(kk))]),
    
    load(['Run2730.6904.',num2str(datasets(kk)),'.mat']);
    
    for i = 1:ncoils
        rd1 = raws2(17,:,:,i);
        rd1 = squeeze(rd1);
        
        rd2 = flipud(rd1(1:floor((size(rd1,1) + 1)/2),:));
        
        rd3 = (rd1(floor((size(rd1,1) + 1)/2 + 1):size(rd1,1),:));
        rd = [rd2; rd3];
        rd = rd.*mask;
        
        rdc1(:,:,i) = rd;        
    end
    
    DATAcomp = rdc1.*compmatrix;
    scale_fctr = norm(DATAcomp(:))/sqrt(ncoils)/20;
    rdc1 = rdc1./scale_fctr;
%     DATAcomp = DATAcomp/scale_fctr;
    
    disp('performing calibration for SPIRiT')
    kCalib = crop(rdc1,[CalibSize, ncoils]);
    kernel = zeros([kSize, ncoils, ncoils]);
    
    [AtA,] = corrMatrix(kCalib,kSize);
    for n=1:ncoils
        kernel(:,:,:,n) = calibrate(AtA,kSize,ncoils,n,CalibTyk);
    end
    GOP = SPIRiT(kernel, 'conv',N);
    
    [res_pocs] = pocsSPIRiT(rdc1,GOP,nIter,rdc1,wavWeight,0);
    
    im_pocsspirit = ifft2c(res_pocs);
    im_pocsspirit_sqr = sos(im_pocsspirit);
    
    subplot(figr,figc,kk);
    imagesc(im_pocsspirit_sqr), title(['Image corresponding to dataset number n = ',num2str(datasets(kk))])
    colormap(gray)
end

