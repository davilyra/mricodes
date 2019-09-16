function res = spirPI(data,weights,traj,imagesize,kernelsize,accel,CalibTyk,lambda)
% res = spirPI(data,weights,traj,imagesize[,kernelsize,accel,CalibTyk,lambda])
%
% this function computes the SPIRiT reconstruction of a MRI data
% Inputs
% - data: parallel MRI data
% - weights: weights of the acquired data
% - traj: sampled trajectory
% - imagesize: the size of the output image
% - kernelsize: which is the size of the kernel
% - accel: acceleration factor that will be used in undersampling data
% - CalibTyk: Tykhonov regularization
% - lambda: ratio between data and calibration consistency.
%           1 recommended when density compensation is used.
%
% This code uses Lustig's SPIRiT toolbox as a base and performs the image
% reconstruction according to his papper
%
% written by Davi Marco Lyra-Leite (davi@ieee.org)
% October 18th 2011

if nargin < 5
    kernelsize = [5,5];         % Standard kernel's size
end

if nargin < 6
    accel = 2;                  % Acceleration factor
end

if nargin < 7
    CalibTyk = 0.0125;          % Tikhonov regularization for calibration 0.01-0.05 recommended
end
if nargin < 8
    lambda = 1;                 % Ratio between data and calibration consistency. 
end

%% Initial Parameters
kSize = kernelsize;               % SPIRiT Kernel size
CalibSize = [30,30];              % size of the calibration region
N = imagesize;                    % size of the target image
nIterCG = 50;                     % number of reconstruction iterations
nCoil = size(data,3);             % number of coils


%% perform SVD coil compression
% D = reshape(data,size(data,1)*size(data,2),size(data,3));
% [~,S,V] = svd(D,'econ');
% nCoil = find(diag(S)/S(1)>0.25, 1, 'last' );
% data = reshape(D*V(:,1:nCoil),size(data,1),size(data,2),nCoil);


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
GOP = SPIRiT(kernel, 'image',N); % init the SPIRiT Operator

% Undersample the data and prepare new NUFFT operator for it
idx = (1:accel:size(traj,2));
k_u = traj(:,idx);
w_u = weights(:,idx);  % use w_u = w(:,idx)*0+1; if you don't want to use density weighting
                       % this may improve noise, but will converge slower. use
                       % larger lambda.
kData_u = data(:,idx,:);

GFFT_u = NUFFT(k_u,w_u, [0,0], N);

im_dc = GFFT_u'*(kData_u.*repmat(sqrt(w_u),[1,1,nCoil]))*accel;
res = im_dc;
res = cgNUSPIRiT(kData_u.*repmat(sqrt(w_u),[1,1,nCoil]),res,GFFT_u,GOP,nIterCG,lambda);
