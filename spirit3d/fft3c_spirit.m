function res = fft3c_spirit(x)
% (c) Michael Lustig 2007
%
% Modified by Taehoon Shin Mar 2011

res = zeros(size(x));
fctr = size(x,1)*size(x,2)*size(x,3);
for n=1:size(x,4)
	%res(:,:,n) = 1/sqrt(fctr)*fftshift(fft2(ifftshift(x(:,:,n))));
	res(:,:,:,n) = 1/sqrt(fctr)*fftshift(fftn(ifftshift( squeeze(x(:,:,:,n)))));
end


