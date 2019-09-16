function [GOP,lambda,nIterCG] = spiritcalibration(kmatrix,coils2recon,undersamplingfactor,iuf,nread,nintl,filename,kxkytraj,kxkyweights,nx,ny)

nVE = size(kmatrix,1);
nphases = size(kmatrix,2);
imagesize = [nx,ny];
intls = 1:nintl;

%disp('starting up the image domain SPIRiT algorithm...')
% Initial Parameters
kSize = [5,5];                       % SPIRiT Kernel size
CalibSize = [30,30];                 % size of the calibration region
CalibTyk = 0.0125;                   % Tykhonov regularization for calibration 0.01-0.05 recommended
lambda = 1;
nIterCG = 50;                     % number of reconstruction iterations

disp(' '),disp(sprintf('starting SPIRiT calibration...'))

kvCalib = nVE/2+1;
ncoils = length(coils2recon);
nphasesCalib = floor(nphases/undersamplingfactor)*undersamplingfactor;
phasesCalib = find(~isnan(kmatrix(kvCalib,1:nphasesCalib)));

% reads data from disk and undersamples it
rkvcpi = nan(nread,1,ncoils,nphasesCalib,nintl);
for p = 1:nphasesCalib,
    k = kmatrix(kvCalib,p);
    if ~isnan(k),
        intlsCalib = (k:iuf:nintl);
        intlszero = setdiff(intls,intlsCalib);
        rawdata = rawloadHD_jfn(filename,[],[],[], 1,coils2recon,p,intlsCalib); %size(rawdata) = [readout kv coil phases intls]
        rkvcpi(:,1,:,p,intlsCalib) = rawdata(1:nread,kvCalib,:,1,:);
        rkvcpi(:,1,:,p,intlszero) = 0;
    else
        rkvcpi(:,1,:,p,:) = 0;
    end;
 end;
kxkyc = undersamplingfactor*permute(mean(rkvcpi,4),[1 5 3 2 4]);

% Gridding and density compensation reconstruction
GFFT1 = NUFFT(kxkytraj,kxkyweights, [0,0] , imagesize);
im = GFFT1'*(kxkyc.*repmat(sqrt(kxkyweights),[1,1,ncoils]));
kData = fft2c(im);
kernel = zeros([kSize,ncoils,ncoils]);

% %%%%%%%%%%%
% subplot(221),imagesc(abs(im(:,:,1))),colorbar,axis image
% subplot(222),imagesc(abs(im(:,:,2))),colorbar,axis image
% subplot(223),imagesc(abs(im(:,:,3))),colorbar,axis image
% subplot(224),imagesc(abs(im(:,:,4))),colorbar,axis image
% colormap gray

if ncoils>1,
    kCalib = crop(kData,[CalibSize,ncoils]); % crop center k-space for calibration
else,
    kCalib = crop(kData,CalibSize);
end;

%prepare calibration matrix for calibration ( Phil Beatty calls this
%correlation values matrix. See thesis by Phil Beatty.)
AtA = corrMatrix(kCalib,kSize);
for n=1:ncoils
    kernel(:,:,:,n) = calibrate(AtA,kSize,ncoils,n,CalibTyk);
end
GOP = SPIRiT(kernel,'image',imagesize); % init the SPIRiT Operator