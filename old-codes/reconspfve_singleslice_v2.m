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
% Changes by Davi Marco Lyra Leite
% October 19, 2011

clear,clc,close all

slice = 4; %prescribe the slice number here (1 to 5)

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
accel = 1;                        % Acceleration factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reconstructs the time-velocity distributions from the selected slice
disp(' '),disp(sprintf('reconstructing data from slice %d...',slice))
filename = sprintf('%sslice%d.7',rawpath,slice); %filename of the current slice
xykvt = zeros(pixels,pixels,nVE,nphases); % dimensions are: x-y-kv-time
xykvt1 = zeros(pixels,pixels,nVE,nphases); % dimensions are: x-y-kv-time
res = zeros(pixels,pixels,nVE,nphases);
for p = 1:nphases, %loops through the cardiac cycle
    xykvc = zeros(pixels,pixels,nVE,ncoils); % kv data in a each voxel for each coil
    for coil = 1:ncoils, %loops through each coil element
        rawdata = rawloadHD_jfn(filename,[],[],[], 1,coil,p,1:nintl); %reads one cardiac phase from each coil
        for v = 1:nVE, %loops through the kv levels
            kxkydata = permute(rawdata(:,v,1,1,:),[1 5 2 4 3]); %rearranges data from each kv level
            weighteddata = kxkydata(:).*kxkyweights(:); %weights the kspace (kx-ky) data based on sampling density
            m = nufft_adj(weighteddata,st)/pixels; %non-cartesian inverse Fourier transform along kx-ky
            xykvc(:,:,v,coil) = m.'; %flips the image and stores it
                        
        end;
    end;
    
    %% SPIRiT implementation
    data = xykvc;
	data = permute(data(:,:,nVE/2+1,:),[1 2 4 3]);
	
    % perform SVD coil compression
    D = reshape(data,size(data,1)*size(data,2),size(data,3));
    [~,S,V] = svd(D,'econ');
    nCoil = find(diag(S)/S(1)>0.25, 1, 'last' );
    data = reshape(D*V(:,1:nCoil),size(data,1),size(data,2),nCoil);

    % Calibration
    % Gridding and density compensation reconstruction
    GFFT1 = NUFFT(kxkytraj,kxkyweights, [0,0] , N);
    im = GFFT1'*(data.*repmat(sqrt(kxkyweights),[1,1,nCoil]));
    kData = fft2c(im);

    kernel = zeros([kSize,ncoils,ncoils]);
    kCalib = crop(data,[CalibSize,ncoils]); % crop center k-space for calibration
            
    [AtA,] = corrMatrix(kCalib,kSize);
    for n=1:ncoils
        kernel(:,:,:,n) = calibrate(AtA,kSize,ncoils,n,CalibTyk);
    end
    GOP = SPIRiT(kernel, 'fft',N); % init the SPIRiT Operator
        
    idx = (1:accel:size(om,2));
    k_u = om(:,idx);
    w_u = kxkyweights(:,idx);  % use w_u = w(:,idx)*0+1; if you don't want to use density weighting
                               % this may improve noise, but will converge slower. use
                               % larger lambda.
	kData_u = data(:,idx,:);

    GFFT_u = NUFFT(k_u,w_u, [0,0], N);

    im_dc = GFFT_u'*(kData_u.*repmat(sqrt(w_u),[1,1,nCoil]))*accel;
    res = im_dc;
    res = cgNUSPIRiT(kData_u.*repmat(sqrt(w_u),[1,1,nCoil]),res,GFFT_u,GOP,nIterCG,lambda);
    [res] = cgSPIRiT(data,GOP, nIterCG, lambda, data);
    
    % returning to xykvc matrix
    a = xykvc(1,1,:,1);
    a = a(:);
    r1 = res(:,1,1);
    r2 = res(1,:,1); r2 = r2(:);
    r3 = res(1,1,:); r3 = r3(:);
    xykvc1 = zeros(size(xykvc));
    xykvc1(:,1,1,1) = r1;
    xykvc1(1,:,1,1) = r2;
    xykvc1(1,1,:,1) = a;
    xykvc1(1,1,1,:) = r3;
        
    clear a r1 r2 r3
    
%     xykvt1(:,:,:,p) = sos(xykvc1,4);
    xykvt(:,:,:,p) = combine4channels(xykvc(:,:,:,1),xykvc(:,:,:,2),xykvc(:,:,:,3),xykvc(:,:,:,4)); %combines the data from the 4 coils
    xykvt1(:,:,:,p) = combine4channels(xykvc1(:,:,:,1),xykvc1(:,:,:,2),xykvc1(:,:,:,3),xykvc1(:,:,:,4)); %combines the data from the 4 coils
    disp(sprintf('%3d/%d',p,nphases)) %display current cardiac phase
end;

xyvt = fftshift(ifft(fftshift(xykvt,3),[],3),3); %the time velocity distsribution is calculated by inverse fourier transform along kv
xyvt1 = fftshift(ifft(fftshift(xykvt1,3),[],3),3); %the time velocity distsribution is calculated by inverse fourier transform along kv

figure,
disp(' ')
disp('Plotting data. Close the figure or press Ctrl-C to stop.')

while(1),
    %plots the magnitude image (central kv, first cardiac phase)
    subplot(211),
    imshow(abs(xykvt1(:,:,nVE/2+1,1)),[ ])
    set(gca,'YDir','normal')
    title(sprintf('slice %d',slice))

    %user will click on a pixel of the image
    xlabel('CLICK ON A BLOOD VESSEL / CLOSE FIGURE TO STOP')
    [Y,X] = ginput(1);

    %takes the velocity distribution from the prescribed pixel
    vt = xyvt1(round(X),round(Y),:,:);
    vt = permute(vt,[3 4 1 2]);

    %plots the time-velocity distribution
    taxis = (0:(nphases-1))*optr/1000;
    vaxis = (((0:(nVE-1))/nVE)-0.5)*maxveloc*2;
    subplot(212),imagesc(taxis,vaxis,abs(vt))
    colormap(gray)
    set(gca,'YDir','normal')
    xlabel('time (ms)')
    ylabel('velocity (cm/s)')

end
