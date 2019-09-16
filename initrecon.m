function [GFFT,wroot,intls2read] = initrecon(ncoils,iuf,kxkytraj,kxkyweights,nx,ny)
% function [GFFT,wroot,intls2read] = initrecon(ncoils,iuf,kxkytraj,kxkyweights,nx,ny)
% This function was written to be used in the reconstruction of spiral
% datasets (especially variable density spirals). It is based on the NUFFT
% algorithm as incorporated in the SPIRiT toolbox, by Miki Lustig.
%
% Inputs:
%   - ncoils: # of coils
%   - iuf: undersampling factor
%   - kxkytraj: k-space trajectory
%   - kxkyweights: k-space weights
%   - [nx,ny]: output image size in the x and y directions, respectively
%
% Outputs:
%   - GFFT: nuFFT object to be used in the image recon
%   - wroot: tensor that contains the squareroot of the k-space weights
%   - intls2read: which spiral interleaves should be read over the original
%   ones
%
% Joao L. A. Carvalho, Aug 2012
%
% Edited by Davi M. Lyra-Leite, Apr 22, 2013.

imagesize = [nx,ny];
nintl = size(kxkytraj,2);

intlsperphase = nintl/iuf;
intls2read = zeros(intlsperphase,iuf);

for k = 1:iuf,
    intls2read(:,k) = (k:iuf:nintl)';
end;

wroot = repmat(sqrt(kxkyweights),[1,1,ncoils]);
GFFT = NUFFT(kxkytraj,kxkyweights,[0,0],imagesize);
