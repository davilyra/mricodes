function [kernel,rawkernel] = calibrate3d(AtA, kSize, nCoil, coil, lambda, sampling)
% function [kernel,rawkernel] = calibrate3d(AtA, kSize, nCoil, coil, lambda, sampling)
% Inputs:
%   - AtA: correlation values matrix, obtained using 
%       [AtA,A] = corrMatrix3d(kCalib, kSize), where kCalib is the
%       auto-calibration set and kSize is the kernel size.
%   - kSize: SPIRiT3D kernel size, an interesting value is: [7,7,3].
%   - nCoil: total number of coils
%   - coil: coil to be calibrated
%   - lambda: Tikhonov regularization parameter (recommended value <= 0.02)
%   - sampling: sampling pattern over the auto-calibration data, it should
%   already be in a Cartesian trajectory
%
%  Outputs:
%   - kernel: calibration kernel (or calibration matrix)
%
% (c) Michael Lustig 2007
%
% Modified by Taehoon Shin Mar 2011
%
% Edited by Davi M. Lyra-Leite, Apr 2013

if nargin < 6
	sampling = ones([kSize,nCoil]);
end

dummyK = zeros(kSize(1),kSize(2),kSize(3), nCoil); 
dummyK((end+1)/2,(end+1)/2,(end+1)/2,coil) = 1;

idxY = find(dummyK);
sampling(idxY) = 0;
idxA = find(sampling);

Aty = AtA(:,idxY); Aty = Aty(idxA);
AtA = AtA(idxA,:); AtA =  AtA(:,idxA);

% kernel = sampling*0;
kernel = zeros(size(sampling));

lambda = (norm(AtA,'fro')/size(AtA,1))*lambda;

% rawkernel = inv(AtA + eye(size(AtA))*lambda)*Aty;
rawkernel = (AtA + eye(size(AtA))*lambda)\Aty;
kernel(idxA) = rawkernel; 
