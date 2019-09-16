function circmask = mask(pixels,radius)

if nargin <2,
    radius = pixels/2;
end;

x = ((-pixels/2):(pixels/2-1))/pixels;
[x y]=meshgrid(x,x);
r = sqrt(x.^2+y.^2);
circmask = r<(radius/pixels);