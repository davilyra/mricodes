function m = idrft1d(d,w,k,npts)
% USAGE: m = idrft1d(d,w,k,npts);
%
% INPUTS: 
% d - vector with k-space samples (raw data)
% k - vector with k-space coordinates, normalized s.t. -0.5<k<0.5
% w - vector with corresponding density-compensation weights (normalized s.t. weight for k=0 is 1)
% npts - number of samples in output vector (must be even)
%
% OUTPUTS:
% m - reconstructed vector
%
% Joao L A Carvalho 2011-01-13 joaoluiz@pgea.unb.br

%density precompensation (weighting)
dw = d.*w;

%total number of k-space samples
N = length(dw);

%clear m(x) before starting
m = zeros(npts,1);

%for each position , compute the Fourier summation of
%the weighted raw data samples
for x=(-npts/2):(npts/2-1),
    m(x+1+npts/2) = sum(dw.*exp(j*2*pi*k*x));
end;