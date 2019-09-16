function range = dispimageD( im, range, nshades );
% range = displayImage (mtx, range, nshades )
% 
% Display a MatLab matrix in the current figure, as a grayscale image.
% 
% RANGE (optional) is a 2-vector specifying the values that map to
% black and white, respectively. Passing a value of 'AUTO1' (default)
% auto-scales to fill range.  'AUTO2' puts image mean at 50% gray.
% 
% NSHADES (optional) specifies the number of gray shades, and defaults
% to the size of the current colormap.
%
% As a modification of dispimage EPS, 6/96 to show the image in its real
% scale in order to not to lose information by MATLAB's distoritions.
%
% DMLL, 01/2011.

if (exist('range') ~= 1)
  range = 'auto1';
end

if (exist('nshades') ~= 1)
  nshades = size(colormap,1);
end

if (nshades < 2)
  nshades = 2;
end

if (strcmp(range,'auto1') || strcmp(range,'auto'))
  if isreal(im)
    mn = min(min(im));
    mx = max(max(im));
  else
    mn = min(min(min(real(im))),min(min(imag(im))));
    mx = max(max(max(real(im))),max(max(imag(im))));
  end
  range = [mn,mx];
elseif strcmp(range,'auto2')
  if isreal(im)
    mn = min(min(im));
    mx = max(max(im));
    av = mean2(im);
  else
    mn = min(min(min(real(im),imag(im))));
    mx = max(max(max(real(im),imag(im))));
    av = (mean2(real(im)) + mean2(imag(im)))/2;
  end
  mx = max(mx-av,av-mn);
  range = [av-mx,av+mx];
end    

if ((range(2) - range(1)) <= eps)
  range(1) = range(1) - 0.5;
  range(2) = range(2) + 0.5;
end

colormap(gray(nshades));

if isreal(im)
  factor=1;
else
  factor = 1+sqrt(-1);
end

%% same as (im-mn)/(mx-mn) + 1:
d_im = im * ((nshades-1) / (range(2)-range(1))) + ...
    factor*(1.5 - (range(1)*(nshades-1))/(range(2)-range(1)));

ii = d_im;

ii = ii/(max(ii(:)));
% ii = ii/(sqrt(sum(abs(ii(:).^2))));
imshow(ii,'InitialMagnification',100);
