function im_ss = sumsquare( im, ind_coil, w)
%  im_ss = sumsquare( im, ind_coil, w)
%
%   Taehoon Shin (shinage@gmail.com) Mar 2008


%ncoils = size( im,1 );

if (nargin < 2 )
    ind_coil = 1:size(im,1);
end

if (nargin < 3 )
    w = ones(1, length(ind_coil));
end



imsize = size( im);
im_ss = zeros( imsize(2:end));


for cc = 1: length(ind_coil);
    im_ss = im_ss + w(cc)*reshape( ( abs( im( ind_coil(cc),:)) ).^2, imsize(2:end) );
end
im_ss = sqrt( im_ss );
