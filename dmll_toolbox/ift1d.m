function m=ift1d(M)
% USAGE: m=ift1d(M)
m = length(M)*fftshift(ifft(fftshift(M)));