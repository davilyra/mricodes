function res = spirPI3(data,rawdata2recon,weights,kweights2recon,traj,ktraj2recon,imagesize,kernelsize,accel)
% res = spirPI(data,weights,traj,imagesize[,kernelsize,accel])
%
% this function computes the SPIRiT reconstruction of a MRI data
% Inputs
% - data: parallel MRI data
% - deights: weights of the acquired data
% - traj: sampled trajectory
% - imagesize: the size of the output image
% - kernelsze: which is the size of the kernel
% - accel: acceleration factor that will be used in undersampling data
%
% This code uses Lustig's SPIRiT toolbox as a base and performs the image
% reconstruction according to his papper
%
% written by Davi Marco Lyra-Leite (davi@ieee.org)
% October 18th 2011

if nargin < 8
    kernelsize = [5,5];            % Standard kernel's size
end

if nargin < 9
    accel = 2;                      % Acceleration factor
end


%% Initial Parameters
kSize = kernelsize;               % SPIRiT Kernel size
CalibSize = [30,30];              % size of the calibration region
N = imagesize;                    % size of the target image
nIterCG = 50;                     % number of reconstruction iterations
CalibTyk = 0.0125;                % Tykhonov regularization for calibration 0.01-0.05 recommended
lambda = 1;                       % Ratio between data and calibration consistency. 1 recommended when density compensation is used.
nCoil = size(data,3);


%% Calibration
% Gridding and density compensation reconstruction
GFFT1 = NUFFT(traj,weights, [0,0] , N);
im = GFFT1'*(data.*repmat(sqrt(weights),[1,1,nCoil]));
kData = fft2c(im);

kernel = zeros([kSize,nCoil,nCoil]);
kCalib = crop(kData,[CalibSize,nCoil]); % crop center k-space for calibration

%prepare calibration matrix for calibration ( Phil Beatty calls this
%correlation values matrix. See thesis by Phil Beatty.)
[AtA,] = corrMatrix(kCalib,kSize);
for n=1:nCoil
	kernel(:,:,:,n) = calibrate(AtA,kSize,nCoil,n,CalibTyk);
end
GOP = SPIRiT(kernel,'image',N); % init the SPIRiT Operator

% % Undersample the data and prepare new NUFFT operator for it
% idx = (1:accel:size(traj,2));
% k_u = traj(:,idx);
% w_u = weights(:,idx);  % use w_u = w(:,idx)*0+1; if you don't want to use density weighting
%                        % this may improve noise, but will converge slower. use
%                        % larger lambda.
% kData_u = data(:,idx,:);

GFFT_u = NUFFT(ktraj2recon,kweights2recon, [0,0], N);

im_dc = GFFT_u'*(rawdata2recon.*repmat(sqrt(kweights2recon),[1,1,nCoil]))*accel;
res = im_dc;
res = cgNUSPIRiT(rawdata2recon.*repmat(sqrt(kweights2recon),[1,1,nCoil]),res,GFFT_u,GOP,nIterCG,lambda);

nIter = 5;
res = pocs((res), GOP, (res), imagesize, nIter, 0);


function X = pocs(data, GOP, x0, imagesize, nIter, show)
% X = spiralpocs(data, GOP, x0, imagesize, nIter, show)
%
% Implementation of the non-Cartesian, POCS l1-SPIRiT reconstruction
%
% Input:
%	- data - Reconstructed data. Make sure that empty entries are zero or else they will not be filled.
%	- GOP: SPIRiT object
%   - x0: initial guess
%	- imagesize: which is the final image size
%	- nIter -	Maximum number of iterations
%	- show: 1 to display progress (slower)
%
% Outputs:
%	- X: reconstructed images for each coil
%
% Based on the original code by Michael Lustig, from 2007
%
% This version by Davi M. Lyra-Leite <davi@ieee.org>, Apr 23, 2013

%==========================================================================
% Initial parameters
wavWeight = 1.5e-2;  % Wavelet soft-thresholding regularization in the reconstruction (SPIRiT only)


%==========================================================================
% find the closest diadic size for the images
sx = imagesize(1);
sy = imagesize(2);
nc = size(data,3);
ssx = 2^ceil(log2(sx)); 
ssy = 2^ceil(log2(sy));
ss = max(ssx, ssy);
W = Wavelet('Daubechies',4,4);
% W = Wavelet('Haar',2,3);

x = x0;
immask = repmat(ones(sx,sy),[1,1,size(data,3)]);

for n=1:nIter
    x = (GOP*x + x).*immask + data;
    
    % apply wavelet thresholding
    X = (x); % goto image domain
    X = zpad(X,ss,ss,nc); % zpad to the closest diadic
    X = W*(X); % apply wavelet
    X = softThresh(X,wavWeight); % threshold ( joint sparsity)
    X = W'*(X); % get back the image
    X = crop(X,sx,sy,nc); % return to the original size
    xx = (X); % go back to k-space
    x = xx.*immask + data; % fix the data
    
    
    if show
        X = ifft2c(x);
        Xsqr = sqrt(sum(abs(X).^2,3));
        figure(show), imshow(Xsqr,[],'InitialMagnification',400); drawnow
    end
end
X = (x); % goto image domain



function x = softThresh(y,t)
    % apply joint sparsity soft-thresholding 
    absy = sqrt(sum(abs(y).^2,3));
    unity = y./(repmat(absy,[1,1,size(y,3)])+eps);

    res = absy - t;
    res = (res + abs(res))/2;
x = unity.*repmat(res,[1,1,size(y,3)]);

