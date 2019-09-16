addpath('/home/mri-vh/Davi/SPIRiT_v0.1')
addpath('/home/mri-vh/Davi/SPIRiT_v0.1/utils')

clear all;
% close all;
clc;

datasets = 12;
% datasets = 10:18;
dts = length(datasets);

% data and recon information
ncoils = 10;        % number of coils
N = [128, 256];     % imagesize


%% Sampling
% pseudorandom sampling
pdf = genPDF(N, 15, 4/8, 2, 0, 0);
mask = genSampling(pdf,500,2);


%% SPIRiT parameters
kSize = [7,7];      % SPIRiT kernel size
nIter = 05;         % number of iteration; phantom requires twice as much as the brain
CalibTyk = 0.1;     % Tykhonov regularization in the calibration
wavWeight = 0.15;   % Wavelet soft-thresholding regularization in the reconstruction


%% CS parameters
TVWeight = 0.005;
xfmWeight = 0.01;
Itnlim = 25;

% Fourier sampling operator
FT = p2DFT(mask, N, 1, 2);
XFM = Wavelet('Daubechies',4,4);
% XFM = 1;

param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;


%% SPIRiT Calibration
[CalibSize, dcomp] = getCalibSize(mask);  % get size of calibration area from mask
compmatrix = repmat(dcomp,[1,1,ncoils]);  % compensation matrix


%% Initializing reconstruction matrices
rdc0 = zeros(N(1), N(2), ncoils); % fully-sampled data matrix
rdc1 = zeros(N(1), N(2), ncoils); % SPIRiT and CS data matrix
rdc2 = zeros(N(1), N(2), ncoils); % CS recon matrix

% figure; clf;
% figr = floor(sqrt(dts));
% figc = ceil(dts/figr);


%% Recon
for kk = 1:dts
    disp(['reconstructing dataset number ',num2str(datasets(kk))]),
    load(['Run2730.6904.',num2str(datasets(kk)),'.mat']);
    
    for i = 1:ncoils
        % Organizing the data
        rd1 = raws2(17,:,:,i);
        rd1 = squeeze(rd1);
        
        rd2 = flipud(rd1(1:floor((size(rd1,1) + 1)/2),:));
        
        rd3 = (rd1(floor((size(rd1,1) + 1)/2 + 1):size(rd1,1),:));
        rd = [rd2; rd3];
        rdc0(:,:,i) = rd; % Fully-sampled data
        
        rd = rd.*mask; % Sampling
        rdc1(:,:,i) = rd; % SPIRiT data
    end
    
    DATAcomp = rdc1.*compmatrix;
    scale_fctr = norm(DATAcomp(:))/sqrt(ncoils)/20;
    rdc0 = rdc0./scale_fctr; % FS Recon
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
    im_pocsspirit_sqr = sos(im_pocsspirit);
    
    
    % CS Recon
    disp('performing CS reconstruction')
    for i = 1:ncoils
        param.data = rdc1(:,:,i);
        im_dc = FT'*(rdc1(:,:,i));
        res = XFM*im_dc;
        
        for n=1:5
            res = fnlCg(res,param);
        end
        
        rdc2(:,:,i) = XFM'*res;
    end
    rdc3 = sos(rdc2);
    rdc = sos(ifft2c(rdc0));
    
    disp('Done!!')
    
    figure;
    subplot(2,3,1),imagesc(abs(rdc)), title('Original sampled data'); colormap(gray)
    subplot(2,3,2),imagesc(abs(im_pocsspirit_sqr)), title('SPIRiT Recon'); colormap(gray)
    subplot(2,3,3),imagesc(abs(rdc3)), title('CS Recon'); colormap(gray)
    subplot(2,3,4),imagesc(abs(sos(ifft2c(rdc1)))), title('SoS'); colormap(gray)
    subplot(2,3,5),imagesc(abs(abs(rdc) - abs(im_pocsspirit_sqr))), title('SPIRiT error'); colormap(gray)
    subplot(2,3,6),imagesc(abs(abs(rdc) - abs(rdc3))), title('CS error'); colormap(gray)
    
%     subplot(figr,figc,kk);
%     imagesc(im_pocsspirit_sqr), title(['Image corresponding to dataset number n = ',num2str(datasets(kk))])
%     colormap(gray)
end

