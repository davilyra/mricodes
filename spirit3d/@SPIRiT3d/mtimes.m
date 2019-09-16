function res = mtimes(a,x)
% This method performs the SPIRiT operator

kernel = a.kernel;  % kernel-[7,7,3,4,4]
nCoils = size(kernel,5);
kSize = size(kernel); kSize = kSize(1:3);
method = a.method;
KERNEL = a.KERNEL;

res = zeros(size(x));
if a.adjoint
    for n=1:nCoils
        tmpk = squeeze(conj(KERNEL(:,:,:,n,:)));
        res(:,:,:,n) = sum(tmpk.*x,4);
    end
    res = res - x;
else
    for n=1:nCoils
        tmpk = KERNEL(:,:,:,:,n);
        res(:,:,:,n) = sum(tmpk.*x,4);
    end
    res = res-x;
end

    
% switch method
% 
% case {'conv'}
% 	res = zeros(size(x));
%     if a.adjoint
%         for n=1:nCoils
%             tmpk = kernel(:,:,:,n);
%             tmpk(floor(kSize(1)/2)+1, floor(kSize(2)/2)+1,n) = ...
% 			tmpk(floor(kSize(1)/2)+1, floor(kSize(2)/2)+1,n) -1;
%             res = res + adjSPIRiT(x(:,:,n),tmpk); 
%         end
%     else
%         for n=1:nCoils
%             tmpk = kernel(:,:,:,n);
%             tmpk(floor(kSize(1)/2)+1, floor(kSize(2)/2)+1,n) = ...
% 		    tmpk(floor(kSize(1)/2)+1, floor(kSize(2)/2)+1,n) -1;
%             res(:,:,n) = SPIRiT(x, tmpk);
%         end
% 
%     end
% 
% case {'fft'}
% 	
% 	res = zeros(size(x));
%     if a.adjoint
%         
%         xx = ifft2c(x);
%         for n=1:nCoils
%             tmpk = squeeze(conj(KERNEL(:,:,n,:)));
%             res(:,:,n) = sum(tmpk.*xx,3); 
%         end
%         res = fft2c(res)-x;
%     else
%         xx = ifft2c(x);
%         for n=1:nCoils
%             tmpk = KERNEL(:,:,:,n);
%             res(:,:,n) = sum(tmpk.*xx,3);
%         end
%         res = fft2c(res)-x;
%     end
%     
%  case {'image'}
%     res = zeros(size(x));
%     if a.adjoint
%         for n=1:nCoils
%             tmpk = squeeze(conj(KERNEL(:,:,n,:)));
%             res(:,:,n) = sum(tmpk.*x,3); 
%         end
%         res = res - x;
%     else
%         for n=1:nCoils
%             tmpk = KERNEL(:,:,:,n);
%             res(:,:,n) = sum(tmpk.*x,3); 
%         end
%         res = res-x;
%     end
%         
% end
