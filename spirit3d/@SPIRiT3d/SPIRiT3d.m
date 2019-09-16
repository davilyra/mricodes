function  res = SPIRiT3d(kernel,method, imSize)
%res = SPIRiT3d(kernel,[method, imSize])
%	Implementation of the SPIRiT kernel operator
%
%   Constructor inputs:
%       kernel: [kx,ky,nCoils, nCoils] the spirit 2D convolution kernel.
%               See corrMatrix.m, calibrate.m and the demos on how to
%               generate a kernel.
%       method: implementation type of the operator. There are three
%               'conv'  is a k-space convolution implementation (slowest)
%                       that operates on k-space data and produces k-space data.
%               'fft'   is an fft based implementation of the k-space
%                       convolution through image domain multiplication. 
%                       It operates on k-sapce data and produces k-space data. 
%               'image' The SPIRiT operator is applied to image space data.
%                       This is useful when using image-based non-cartesian
%                       reconstruction. In this case the SPIRiT operator
%                       operates on image data and produces image data. 
%
%       imSize:  Size of the resulting image (only needed for 'fft' and
%                'image' modes
%
%
%   Example:
%
%   See Lustig's example for 2D SPIRiT.
%   
%   See Also:
%               calibrate3d.m, corrMatrix3d.m
%
% (c) Michael Lustig 2007
%
% (c) Taehoon Shin 2011
%
% Modified by Davi M. Lyra-Leite, Apr 2013


KERNEL = zeros([imSize(1), imSize(2), imSize(3), size(kernel,4),size(kernel,5)]);
for n=1:size(kernel,5)
%     KERNEL(:,:,:,n) = ifft2c(   zpad(  kernel(end:-1:1,end:-1:1,:,n)*sqrt(imSize(1)*imSize(2)), imSize(1), imSize(2), size(kernel,3)  )    );
    KERNEL(:,:,:,:,n) = ifft3c_spirit( zpad( kernel(end:-1:1,end:-1:1,end:-1:1,:,n)*sqrt(imSize(1)*imSize(2)*imSize(3)), imSize(1), imSize(2), imSize(3), size(kernel,4) )  );
end

res.kernel = kernel;
res.adjoint = 0;
res.KERNEL = KERNEL;
res.method = method;
res.imSize = imSize;
res = class(res,'SPIRiT3d');
