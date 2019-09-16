function dispimage(m,range)

if nargin==2,
    imagesc(m,range),
else
    imagesc(m),
end;
    
axis image off,colormap(gray)