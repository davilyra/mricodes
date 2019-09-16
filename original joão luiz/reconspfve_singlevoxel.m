% RECONSTRUCTION OF THE TIME-VELOCITY DISTRIBUTION FROM
% A SINGLE VOXEL IN A 3D VOLUME, WHICH IS MANUALLY
% PRESCRIBED BY THE USER.
%
% MEMORY REQUIREMENTS ARE NOT SIGNIFICANT.
%
% Written by Joao L. A. Carvalho <joaoluiz@gmail.com>
% Department of Electrical Engineering
% University of Brasilia, Brazil
%
% July 14, 2008

clear,clc,close all

%this sets up fesslers code for the non-Cartesian FFT (nuFFT)
% cd fessler
% setup
% cd ..

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

% This loop reads one image from each slice.
% The image corresponds to the central kv level of the first cardiac phase.
xyzc = zeros(pixels,pixels,nslices,ncoils); %dimensions are: x-y-z-coil
for slice = 1:nslices,    
    filename = sprintf('%sslice%d.7',rawpath,slice); %filename for each slice
    for coil = 1:ncoils,
        rawdata = rawloadHD_jfn(filename,[],[],[], 1,coil,1,1:nintl); %reads the first cardiac phase of each coil
        kxkydata = permute(rawdata(:,nVE/2+1,1,1,:),[1 5 2 4 3]); %rearranges data from the central kv level (nVE/2+1)
        weighteddata = kxkydata(:).*kxkyweights(:); %weights the kspace data based on sampling density
        m = nufft_adj(weighteddata,st)/pixels; %non-cartesian inverse Fourier transform along kx-ky      
        xyzc(:,:,slice,coil) = m.'; %flips the image and stores it in xyzc
    end;   
end;

%combines the data from the 4 coils into a single 3D volume
xyz = combine4channels(xyzc(:,:,:,1),xyzc(:,:,:,2),xyzc(:,:,:,3),xyzc(:,:,:,4));

%rearranges the 3D data into a single 2D image
allslices = zeros(nslices*pixels,pixels);
for slice = 1:nslices,
    x1 = (slice-1)*pixels + 1;
    x2 = x1 - 1 + pixels;
    allslices(x1:x2,:) = xyz(:,:,slice);
end;

%plots the slices
figure,
subplot(121),imshow(abs(allslices),[ ]),
set(gca,'YDir','normal')

%user will click on a slice
title('CLICK ON ONE OF THE SLICES')
[X,Y] = ginput(1);
slice = ceil(Y/pixels);
title(' ')

%plots the selected slice
subplot(222),
imshow(abs(xyz(:,:,slice)),[ ])
set(gca,'YDir','normal')
title(sprintf('slice %d',slice))

%user will click on a pixel of the slice
xlabel('CLICK ON A BLOOD VESSEL')
[Y,X] = ginput(1);

%reconstructs the time-velocity distribution from the selected pixel
xlabel('reconstructing data...'),drawnow
disp(' '),disp('reconstructing data...')
filename = sprintf('%sslice%d.7',rawpath,slice); %filename of the current slice
kvt = zeros(nVE,nphases); % k-t data (kv as a function of time)
for p = 1:nphases, %loops through the cardiac cycle
    kvc = zeros(nVE,ncoils); % kv data in a single voxel (x,y,z) for each coil
    for coil = 1:ncoils, %loops through each coil element
        rawdata = rawloadHD_jfn(filename,[],[],[], 1,coil,p,1:nintl); %reads one cardiac phase from each coil
        for v = 1:nVE, %loops through the kv levels
            kxkydata = permute(rawdata(:,v,1,1,:),[1 5 2 4 3]); %rearranges data from each kv level
            weighteddata = kxkydata(:).*kxkyweights(:); %weights the kspace (kx-ky) data based on sampling density
            m = nufft_adj(weighteddata,st)/pixels; %non-cartesian inverse Fourier transform along kx-ky
            kvc(v,coil) = m(round(Y),round(X)); %takes only the prescribed pixel (note that X and Y are flipped here too)
        end;        
    end;
    kvt(:,p) = combine4channels(kvc(:,1),kvc(:,2),kvc(:,3),kvc(:,4)); %combines the data from the 4 coils
    disp(sprintf('%3d/%d',p,nphases)) %display current cardiac phase
end;
vt = fftshift(ifft(fftshift(kvt,1),[],1),1); %the time velocity distsribution is calculated by inverse fourier transform along kv
xlabel(' ')

%plots the time-velocity distribution
taxis = (0:(nphases-1))*optr/1000;
vaxis = (((0:(nVE-1))/nVE)-0.5)*maxveloc*2;
subplot(224),imagesc(taxis,vaxis,abs(vt))
colormap(gray)
set(gca,'YDir','normal')
xlabel('time (ms)')
ylabel('velocity (cm/s)')



