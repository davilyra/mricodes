function [GFFT,wroot,intls2read]=initrecon(ncoils,iuf,kxkytraj,kxkyweights,nx,ny)

imagesize = [nx,ny];
nintl = size(kxkytraj,2);

intlsperphase = nintl/iuf;
intls2read=zeros(intlsperphase,iuf);
for k = 1:iuf,
    intls2read(:,k) = (k:iuf:nintl)';
end;

wroot = repmat(sqrt(kxkyweights),[1,1,ncoils]);
GFFT = NUFFT(kxkytraj,kxkyweights,[0,0],imagesize);
