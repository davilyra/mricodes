clear,clc,close all

disp('loading data...')
load wholevolume.mat

pixels = size(xyzvt,1);
nslices = size(xyzvt,3);
nVE = size(xyzvt,4);
nphases = size(xyzvt,5);

%rearranges the 3D data into a single 2D image
allslices = zeros(nslices*pixels,pixels);
for slice = 1:nslices,
    x1 = (slice-1)*pixels + 1;
    x2 = x1 - 1 + pixels;
    allslices(x1:x2,:) = xyz(:,:,slice);
end;

disp('Plotting data. Close the figure or press Ctrl-C to stop.')
figure,

while(1),

    %plots the slices
    subplot(121),imshow(allslices,[ ]),
    set(gca,'YDir','normal')

    %user will click on a slice
    title('CLICK ON ONE OF THE SLICES')
    [X,Y] = ginput(1);
    slice = ceil(Y/pixels);
    title(' ')
    
    %plots the selected slice
    subplot(222),
    imshow(xyz(:,:,slice),[ ])
    set(gca,'YDir','normal')
    xlabel(sprintf('slice %d',slice))

    %user will click on a pixel of the slice
    title('CLICK ON A BLOOD VESSEL')
    [Y,X] = ginput(1);
    
    %takes the time-velocity distribution from the prescribed voxel
    vt = xyzvt(round(X),round(Y),slice,:,:);
    vt = permute(vt,[4 5 1 2 3]);

    %plots the time-velocity distribution
    taxis = (0:(nphases-1))*optr/1000;
    vaxis = (((0:(nVE-1))/nVE)-0.5)*maxveloc*2;
    subplot(222),imagesc(taxis,vaxis,abs(vt))
    colormap(gray)
    set(gca,'YDir','normal')
    xlabel('time (ms)')
    ylabel('velocity (cm/s)')

end;