function X = spirPI2(rawdata,kData_u,wroot,w_u,imagesize,undersamplingfactor,CalibTyk,ktraj2recon,show)
% X = spirPI2(rawdata,kData_u,GFFT,weights,w_u,imagesize,undersamplingfactor,CalibTyk)
%
% this function computes the pocsSPIRiT reconstruction of a MRI data
%
% Inputs
% - data: parallel MRI data
% - weights: weights of the acquired data
% - traj: sampled trajectory
% - imagesize: the size of the output image
% - kernelsize: which is the size of the kernel
% - undersamplingfactor: acceleration factor that will be used in undersampling data
% - CalibTyk: Tykhonov regularization
% - lambda: ratio between data and calibration consistency.
%           1 recommended when density compensation is used.
%
% This code uses Lustig's SPIRiT toolbox as a base and performs the image
% reconstruction according to his papper
%
% written by Davi Marco Lyra-Leite (davi@ieee.org)
% October 18th 2011

%% Initial Parameters
kSize = [5,5];                    % SPIRiT Kernel size
CalibSize = [30,30];              % size of the calibration region
nCoil = size(rawdata,3);          % number of coils
nIter = 50;                       % number of iterations


%% Calibration
% Gridding and density compensation reconstruction
w_data = rawdata.*wroot;
im = zeros(imagesize(1),imagesize(1),nCoil);
for i = 1:nCoil
    im(:,:,i) = grid2_deap(w_data(:,:,i),ktraj2recon,imagesize(1));
end

% im = GFFT'*(wroot.*rawdata);
kData = fft2c(im);

kernel = zeros([kSize,nCoil,nCoil]);
kCalib = crop(kData,[CalibSize,nCoil]); % crop center k-space for calibration

% correlation values matrix
[AtA,] = corrMatrix(kCalib,kSize);
for n=1:nCoil
	kernel(:,:,:,n) = calibrate(AtA,kSize,nCoil,n,CalibTyk);
end
GOP = SPIRiT(kernel, 'conv',imagesize); % init the SPIRiT Operator

% Undersample the data and prepare new NUFFT operator for it
% size(kData_u),size(repmat(sqrt(w_u),[1,1,nCoil]))
GFFT_u = NUFFT(ktraj2recon,w_u, [0,0], imagesize);

im_dc = GFFT_u'*(kData_u.*repmat(sqrt(w_u),[1,1,nCoil]))*undersamplingfactor;
res = im_dc;

X = spiralpocs(fft2c(res), GOP, fft2c(res), imagesize, nIter, show);


function X = spiralpocs(data, GOP, x0, imagesize, nIter, show)
% X = spiralpocs(data, GOP, x0, imagesize, nIter, show)
%
% Implementation of the non-Cartesian, POCS l1-SPIRiT reconstruction
%
% Input:
%	- data - Undersampled k-space data. Make sure that empty entries are zero or else they will not be filled.
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
% This version was developed by Davi M. Lyra-Leite <davi@ieee.org>
% Apr 23, 2013

%==========================================================================
% Initial parameters
wavWeight = 0.0015;  % Wavelet soft-thresholding regularization in the reconstruction (SPIRiT only)


%==========================================================================
% find the closest diadic size for the images
sx = imagesize(1);
sy = imagesize(2);
nc = size(data,3);
ssx = 2^ceil(log2(sx)); 
ssy = 2^ceil(log2(sy));
ss = max(ssx, ssy);
W = Wavelet('Daubechies',4,4);
%W = Wavelet('Haar',2,3);

mask = (data == 0);
x = x0;

for n=1:nIter
    x = (x + GOP*x ).*(mask) + data; % Apply (G-I)*x + x
    
    % apply wavelet thresholding
    X = ifft2c(x); % goto image domain
    X= zpad(X,ss,ss,nc); % zpad to the closest diadic
    X = W*(X); % apply wavelet
    X = softThresh(X,wavWeight); % threshold ( joint sparsity)
    X = W'*(X); % get back the image
    X = crop(X,sx,sy,nc); % return to the original size
    xx = fft2c(X); % go back to k-space
    x = xx.*mask + data; % fix the data
    
    if show
        X = ifft2c(x);
        Xsqr = sqrt(sum(abs(X).^2,3));
        figure(show), imshow(Xsqr,[],'InitialMagnification',400);, drawnow
    end
    X = ifft2c(x); % goto image domain
end


function x = softThresh(y,t)
    % apply joint sparsity soft-thresholding 
    absy = sqrt(sum(abs(y).^2,3));
    unity = y./(repmat(absy,[1,1,size(y,3)])+eps);

    res = absy - t;
    res = (res + abs(res))/2;
x = unity.*repmat(res,[1,1,size(y,3)]);
