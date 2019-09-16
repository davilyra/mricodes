% RECONSTRUCTION OF THE TIME-VELOCITY DISTRIBUTIONS FROM
% THE ENTIRE 3D VOLUME.
% THE SPATIAL POSITION IS THEN MANUALLY PRESCRIBED
% BY THE USER FOR PLOTTING. CLOSE THE FIGURE OR PRESS CTRL-C TO STOP.
%
% MATLAB WILL USE 1.7GB OF RAM.
%
% Written by Joao L. A. Carvalho <joaoluiz@gmail.com>
% Department of Electrical Engineering
% University of Brasilia, Brazil
%
% Changes by Davi Marco Lyra-Leite (davi@ieee.org)
% Medical Imaging Research Group
% Department of Electrical Engineering
% University of Brasilia, Brazil
% october 17th 2011

clear,clc,close all

% considering that some important packages are already added to the general
% path of MATLAB (such as Fessler's ir toolbox and Lustig's SPIRiT one), it
% is not necessary to perform their initialization before start the code
% itself

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rawpath = './rawdata/';
nslices = 5;
ncoils = 4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this loads the usercv variables from the rawdata file
filename = sprintf('%sslice1.7',rawpath);
[rawdata,usercv,hdr] = rawloadHD_jfn(filename,[],[],[], 1, [], [], []);
maxveloc = usercv(1); % maximum velocity value (1/2 of velocity field-of-view)
optr = usercv(2); % TR duration (in microseconds)
nphases = usercv(4); %number of cardiac phases (i.e., temporal frames)
nVE = usercv(7); %number of velocity encoding steps
nread = usercv(8); %number of readout samples
nintl = usercv(14); %number of spiral interleaves
spiralid = usercv(10); %this number identifies which spiral trajectory was used
heartrate = usercv(5); %heart rate of the subject during the scan
RRpct = usercv(6); %percentage of cardiac cycle that was covered
flipangle = usercv(15); %flip angle (in degrees)
% vesperbeat = usercv(3);
% vres = usercv(9);
% vegrad = usercv(11);
% realtime = usercv(12);
% nocine = usercv(13);
% ia_rf1 = usercv(16); 
% densityreductionfactor = usercv(17);
% oblique = usercv(18);
% intlsperbeat = usercv(19);

% load spatial parameters
switch(spiralid)
    case 19,
        kfile='recon16cm14mm8intl4vd';
        spatfov = 160;
        spatres = 1.4;         
    otherwise,
        error(['unexpected spiralid value: ',num2str(spiralid)]);
end;
[kxkytraj kxkyweights] = kkread(kfile,nintl,nread);
%figure,plot(kxkytraj);axis equal; %plots kspace trajectory
kxkytraj = kxkytraj / (2*max(abs(kxkytraj(:)))); %normalize to [-0.5 , 0.5]
kxkyweights=kxkyweights/max(kxkyweights(:))/sqrt(2); %normalizes the density values to [0 , sqrt(2)/2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STARTING UP THE NUFFT ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('nufft_init...')
pixels = ceil(spatfov/spatres); %image size
Nd = [pixels,pixels];
Jd = [6,6]; %kernel size (?)
overgridfactor = 2;
om(:,1) = 2*pi*real(kxkytraj(:)); %kspace trajectory
om(:,2) = 2*pi*imag(kxkytraj(:));
st = nufft_init(om, Nd, Jd, overgridfactor*Nd, Nd/2); %initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xykvc = zeros(pixels,pixels,nVE,ncoils); %dimensions are: x-y-kv-coil
xyzvt = zeros(pixels,pixels,nslices,nVE,nphases); %dimensions are: x-y-z-v-t
xyz = zeros(pixels,pixels,nslices); %magnitude images (dimensions are: x-y-z)

disp(' '),disp('reconstructing data...')

for p = 1:nphases, %loops through the cardiac phases
    
    for slice = 1:nslices, %loops through the slices
        filename = sprintf('%sslice%d.7',rawpath,slice); %filename for each slice
        
        for coil = 1:ncoils, %loops through the coil elements
            kxkykv = rawloadHD_jfn(filename,[],[],[], 1,coil,p,1:nintl); %reads one cardiac phase from each coil            
            for v = 1:nVE, %loops through the kv levels                
                kxkydata = permute(kxkykv(:,v,1,1,:),[1 5 2 3 4]); %rearranges data from each kv level
                weighteddata = kxkydata(:).*kxkyweights(:); %weights the kspace data based on sampling density
                xykvc(:,:,v,coil) = nufft_adj(weighteddata,st)/pixels; %non-cartesian inverse Fourier transform along kx-ky
            end;
        end;
            
        xykv = combine4channels(xykvc(:,:,:,1),xykvc(:,:,:,2),xykvc(:,:,:,3),xykvc(:,:,:,4)); %combines the data from the 4 coils
        xyv = fftshift(ifft(fftshift(xykv,3),[],3),3); %the velocity distributions are calculated by inverse fourier transform along kv
        xyzvt(:,:,slice,:,p) = permute(xyv,[1 2 4 3 5]); %rearranges the data and stores the velocity distributions in xyzvt 
        
        if(p==1), %do this only for the first cardiac phase
            xyz(:,:,slice) = abs(xykv(:,:,nVE/2+1)); %takes the central kv level (kv=0)
        end;
            
    end;
    
    disp(sprintf('%3d/%d',p,nphases)) %display current cardiac phase
end;

%flips x and y
xyzvt = permute(xyzvt,[2 1 3 4 5]);
xyz = permute(xyz,[2 1 3]);

%saves the data
disp('saving reconstructed data to disk...')
save wholevolume.mat xyzvt xyz maxveloc optr
disp('Done!')

%rearranges the 3D data into a single 2D image
allslices = zeros(nslices*pixels,pixels);
for slice = 1:nslices,
    x1 = (slice-1)*pixels + 1;
    x2 = x1 - 1 + pixels;
    allslices(x1:x2,:) = xyz(:,:,slice);
end;

figure,
disp(' ')
disp('Plotting data. Close the figure or press Ctrl-C to stop.')

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