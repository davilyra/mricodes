% RECONSTRUCTION OF THE TIME-VELOCITY DISTRIBUTION FROM
% MULTIPLE VOXELS IN A SLICE, WHICH ARE MANUALLY
% PRESCRIBED BY THE USER. CLOSE THE FIGURE OR PRESS CTRL-C TO STOP.
%
% MATLAB WILL USE APPROXIMATELY 900MB OF RAM.
%
% Written by Joao L. A. Carvalho <joaoluiz@gmail.com>
% Department of Electrical Engineering
% University of Brasilia, Brazil
%
% July 14, 2008
%
% Changes by Davi Marco Lyra Leite <davi@ieee.org> e Joao Luiz Carvalho
% November 3rd, 2011

clear,clc,close all

slice = 5; %prescribe the slice number here (1 to 5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rawpath = './rawdata/';
nslices = 5;
ncoils = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this loads the usercv variables from the rawdata file
filename = sprintf('%sslice%d.7',rawpath,slice);
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
% figure,plot(kxkytraj);axis equal; %plots kspace trajectory
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPIRiT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Parameters
kSize = [5,5];                       % SPIRiT Kernel size
CalibSize = [30,30];              % size of the calibration region
N = Nd;                           % size of the target image
nIterCG = 50;                     % number of reconstruction iterations
CalibTyk = 0.02;                  % Tykhonov regularization for calibration 0.01-0.05 recommended
lambda = 1;                       % Ratio between data and calibration consistency. 1 recommended when density compensation is used.
accel = 4;                        % Acceleration factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reconstructs the time-velocity distributions from the selected slice
disp(' '),disp(sprintf('reconstructing data from slice %d...',slice))
filename = sprintf('%sslice%d.7',rawpath,slice); %filename of the current slice
xykvt = zeros(pixels,pixels,1,1); % dimensions are: x-y-kv-time
res = zeros(pixels,pixels,1,1);
p=1;
v=nVE/2+1;
xykvc = zeros(pixels,pixels,1,ncoils); % kv data in a each voxel for each coil
rawdata = rawloadHD_jfn(filename,[],[],[], 1,1:ncoils,p,1:nintl); %reads one cardiac phase from each coil
kxkyckv = permute(rawdata(:,v,:,:,:),[1 5 3 2 4]);

xyc = spirPI(kxkyckv,kxkyweights,kxkytraj,N,kSize,accel);
xykvc(:,:,1,:) = permute(xyc,[1 2 4 3]); % = m.'; %flips the image and stores it
disp(sprintf('kv=%2d/%d p=%3d/%d',v,nVE,p,nphases)) %display current cardiac phase


xykvt(:,:,:,p) = combine4channels(xykvc(:,:,:,1),xykvc(:,:,:,2),xykvc(:,:,:,3),xykvc(:,:,:,4)); %combines the data from the 4 coils
xy = xykvt;

% energia_pi=sum(abs(xy(:).^2)),

%saves the data
cd datapaper
disp('saving reconstructed data to disk...')
save(sprintf('slice_%d_%d_xy.mat',slice,accel), 'xy');
disp('Done!')
cd ..

figure,
imshow(abs(xy),[ ])
set(gca,'YDir','normal')