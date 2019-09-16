function [GOP,lambda,nIterCG] = spiritcalibration2(kmatrix,coils2recon,undersamplingfactor,iuf,nread,nintl,filename,kxkytraj,kxkyweights,nx,ny)

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

vCalib = nVE/2+1;
ncoils = length(coils2recon);
nphasesCalib = floor(nphases/undersamplingfactor)*undersamplingfactor;

% reads data from disk and undersamples it
rkvc_i = zeros(nread,nVE,ncoils,1,nintl);
r_cpi = zeros(nread,1,ncoils,nphasesCalib,nintl);
for p = 1:nphasesCalib,
    rawdata = rawloadHD_jfn(filename,[],[],[], 1,coils2recon,p,intls); %size(rawdata) = [readout kv coil phases intls]
    for v = 1:nVE
        k = kmatrix(v,p);
        if ~isnan(k),
            intlsCalib = k:iuf:nintl;
            rkvc_i(:,v,:,:,intlsCalib) = rawdata(1:nread,v,:,:,intlsCalib);
        end;
    end;
    rvc_i = fftshift(ifft(fftshift(rkvc_i,2),[],2),2);
    r_cpi(:,:,:,p,:) = rvc_i(:,vCalib,:,:,:);
end;
%r_c_i = undersamplingfactor*mean(r_cpi,4);
r_c_i = mean(r_cpi,4);
ric = permute(r_c_i,[1 5 3 2 4]);

% image reconstruction
GFFT = NUFFT(kxkytraj,kxkyweights,[0,0],imagesize);
im = GFFT'*(ric.*repmat(sqrt(kxkyweights),[1,1,ncoils]));

% %%%%%%%%%%%
% subplot(221),imagesc(abs(im(:,:,1))),colorbar,axis image
% subplot(222),imagesc(abs(im(:,:,2))),colorbar,axis image
% subplot(223),imagesc(abs(im(:,:,3))),colorbar,axis image
% subplot(224),imagesc(abs(im(:,:,4))),colorbar,axis image
% colormap gray

kData = fft2c(im);
if ncoils>1,
    kCalib = crop(kData,[CalibSize,ncoils]); % crop center k-space for calibration
else,
    kCalib = crop(kData,CalibSize);
end;

%prepare calibration matrix for calibration ( Phil Beatty calls this
%correlation values matrix. See thesis by Phil Beatty.)
AtA = corrMatrix(kCalib,kSize);
kernel = zeros([kSize,ncoils,ncoils]);
for n=1:ncoils
    kernel(:,:,:,n) = calibrate(AtA,kSize,ncoils,n,CalibTyk);
end
GOP = SPIRiT(kernel,'image',imagesize); % initialize the SPIRiT operator