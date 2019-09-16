function d = idrft2d(m,k)
% USAGE: d = idrft2d(m,k);
%
% INPUTS: 
% k - vector with complex k-space coordinates
% m - Nim x Nim image
%
% OUTPUTS:
% d - vector with k-space samples
%
% Joao L A Carvalho 2010-12-17 joaoluiz@pgea.unb.br

N = length(k); %total number of k-space samples
[nx,ny] = size(m); %image size
d=zeros(size(k)); % initialize k-space data vector

%separates the comples k-space coordinates in kx and ky
kx = real(k);
ky = imag(k);

% d_n = sum_x sum_y m(x,y) exp[-j 2 pi (x k_x + y k_y)]
for x=(-nx/2):(nx/2-1),
    for y=(-ny/2):(ny/2-1),
        d = d+m(x+1+nx/2,y+1+ny/2)*exp(-j*2*pi*(x.*kx+y.*ky));
    end;
end;

%scale factor
d = d/sqrt(nx*ny);