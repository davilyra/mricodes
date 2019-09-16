function [GOP, sampidx_tot,FT_SPR3d,f0] = spiritcalibration3d_davi(kmatrix,coils2recon,undersamplingfactor,iuf,filename,kxkytraj,kxkyweights,imSize,kSize,CalibSize,CalibTyk)
% [GOP, sampidx_tot,FT_SPR3d,f0] = spiritcalibration3d_davi(kmatrix,coils2recon,undersamplingfactor,iuf,filename,kxkytraj,kxkyweights,imSize,[kSize,CalibSize,CalibTyk])
% Performs the spiral SPIRiT 2D calibration
% Inputs:
%   - kmatrix: view ordering scheme in k-space
%   - coils2recon: which coils are wanted to be used to recon the signal
%   - undersamplingfactor: factor that describes how the data is
%   undersampled
%   - iuf: undersampling factor of the interleaves
%   - nread: number of points in the readout
%   - nintl: number of intereaves
%   - filename: the name of the file where data is stored
%   - kxkytraj: trajectory in kx-ky space
%   - kxkyweights: the weights in k-space
%   - imSize of the output image
%   - kSize: SPIRiT3D kernel size
%   - CalibSize: auto-calibration data size
%   - CalibTyk: Tikhonov's regularization factor
% Outputs:
%   - GOP: the SPIRiT initializator
%   - sampidx_tot: Taehoon's sampling matrix
%   - FT_SPR3d: SPIRiT 3D Operator
%   - f0: zero vector with size 1x(#kvs or #kzs)
%
% by Davi Marco Lyra-Leite (davi@ieee.org)
% October 3rd, 2012


if nargin < 11
    CalibTyk = 0.02;
end

if nargin < 10
    CalibSize = [30,30,length(coils2recon)];
end

if nargin < 9
    kSize = [7,7,3];
end

% Using velocity encoding and defining some parameters
nVE = size(kmatrix,1);
nphases = size(kmatrix,2);
nintl = size(kxkytraj,2);
nread = size(kxkytraj,1);
intls = 1:nintl;
ncoils = length(coils2recon);
nphasesCalib = floor(nphases/undersamplingfactor)*undersamplingfactor;
% vCalib = nVE/2+1;

% Starting SPIRiT calibration
disp(' '),disp(sprintf('starting SPIRiT calibration...'))

% reads data from disk and undersamples it
rkvc_i = zeros(nread,nVE,ncoils,1,nintl);
for p = 1:nphasesCalib,
    rawdata = rawloadHD_jfn(filename,[],[],[], 1,coils2recon,p,intls); %size(rawdata) = [nread kv coil phases intls]
    for v = 1:nVE
        k = kmatrix(v,p);
        if ~isnan(k),
            intlsCalib = k:iuf:nintl;
            rkvc_i(:,v,:,:,intlsCalib) = rawdata(1:nread,v,:,:,intlsCalib);
        end
    end
end
% undersampled data
ric = permute(rkvc_i,[1 5 2 3 4]); %size(ric) = [nread nintl kv coil phases]

% f0 = zeros(1,nVE);
f0 = ones(1,nVE);

% sampling matrix according to Taehoon's original code
% matrix size = (number of interleaves effectively used, number of
% elements in kz or kv)
sampidx_tot = ones((nintl/undersamplingfactor),nVE);
for i = 2:(nintl/undersamplingfactor)
    sampidx_tot(i,:) = sampidx_tot(i - 1,:) + undersamplingfactor;
end

for i = 2:undersamplingfactor
    sampidx_tot(:,i:undersamplingfactor:end) = sampidx_tot(:,i:undersamplingfactor:end) + (i - 1);
end


for zz = 1:nVE     
    FT{zz} = NUFFT(kxkytraj,kxkyweights,[0,0],[imSize(1),imSize(2)] );
end

% SPIRiT object (similar to the NUFFT object GFFT1 used in the 2D problem)
FT_SPR3d = SPR3dFULL( FT, [imSize(1),imSize(2),nVE], [nread,nintl,nVE], ncoils, sampidx_tot, f0 );

% calibration target
im_calib = FT_SPR3d'*(ric.*repmat(sqrt(kxkyweights),[1,1,nVE,ncoils]));

kData = fft3c_spirit(im_calib);

if ncoils > 1
	kCalib = crop3d(kData,[CalibSize,ncoils]); % crop center k-space for calibration
else
    kCalib = crop(kData,CalibSize);
end
clear kData;

% prepare calibration matrix for calibration
[AtA,] = corrMatrix3d(kCalib,kSize);

kernel = zeros([kSize,ncoils, ncoils]);
for cc = 1:ncoils
%     disp(sprintf('Calibrating coil %d',cc));
    [kernel_cc ~] = calibrate3d( AtA, kSize, ncoils, cc, CalibTyk);
    kernel(:,:,:,:,cc) = kernel_cc;        
end
disp('Done Calibrating');

% SPIRiT itializator
GOP = SPIRiT3d(kernel, 'image', [imSize(1) imSize(2) nVE] );
