function [GFFT,wroot,intls2read,w_u] = initrecon3(ncoils,iuf,kxkytraj,kxkyweights,nx,ny,nVE)
% function [GFFT,wroot,intls2read,w_u] = initrecon3(ncoils,iuf,kxkytraj,kxkyweights,nx,ny,nVE)
% Inputs:
%   - ncoils: # of coils
%   - iuf: undersampling factor
%   - kxkytraj: k-space trajectory
%   - kxkyweights: k-space weights
%   - [nx,ny]: output image size in the x and y direction, respectivelly
%   - nVE: # of velocity encodings
%
% Outputs:
%   - GFFT:
%   - wroot: tensor containing the square-root of the weights
%   - intls2read: which interleaves are being considered in the
%   retrospectively undersampling
%   - w_u: undersampled weights

imagesize = [nx,ny];
nintl = size(kxkytraj,2);

intlsperphase = nintl/iuf;
intls2read = zeros(intlsperphase,iuf);

for k = 1:iuf,
    intls2read(:,k) = (k:iuf:nintl)';
end;

wroot = repmat(sqrt(kxkyweights),[1,1,ncoils,nVE]);
% GFFT = NUFFT(kxkytraj,kxkyweights,[0,0],imagesize);

idx = (1:iuf:size(kxkytraj,2));
k_u = zeros(size(kxkytraj));
k_u(:,idx) = kxkytraj(:,idx);
w_u = zeros(size(kxkyweights));
w_u(:,idx) = kxkyweights(:,idx);  % use w_u = w(:,idx)*0+1; if you don't want to use density weighting
                       % this may improve noise, but will converge slower. use
                       % larger lambda.

for zz = 1:nVE   
    GFFT{zz} = NUFFT(k_u,w_u,[0,0],[ imagesize(1),imagesize(2)] );
end