function m = idrft(d,w,k,Nim)
% USAGE: m = idrft(d,w,k,Nim);
%
% INPUTS: 
% d - vector with k-space samples
% w - vector with corresponding weighting
% k - vector with complex k-space coordinates
% Nim - image size in pixels
%
% OUTPUTS:
% m - Nim x Nim reconstructed image
%
% Joao L A Carvalho 2010-12-17 joaoluiz@pgea.unb.br

%image size
nx = Nim;
ny = Nim;

%density precompensation (weighting)
dw = d.*w;

%total number of k-space samples
N = length(dw);

%separates the comples k-space coordinates in kx and ky
kx = real(k);
ky = imag(k);

%clear m(x,y) before starting
m = zeros(nx,nx);

%for each pixel (x,y), compute the Fourier summation of
%the weighted raw data samples
for x=(-nx/2):(nx/2-1),
    for y=(-ny/2):(ny/2-1),
        m(x+1+nx/2,y+1+ny/2) = sum(dw.*exp(j*2*pi*(x.*kx+y.*ky)));
    end;
end;

%scale factor
m = m/sqrt(nx*ny);
