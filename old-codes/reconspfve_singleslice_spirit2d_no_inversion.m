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
% October 31, 2011

%% starting the code by cleaning memory and closing opened figures
clear,clc,close all

%% setting code parameters
savedata2disk = 0; % [0] no, [1] yes
rawpath = './rawdata/';
slices2recon = [1]; %examples: [4] or [1 2 3 4 5]
coils2recon = [1 2 3 4]; %[1 2 3 4]
ncoils = 4;

slice = 4; %prescribe the slice number here (1 to 5)

%==========================================================================
% INITIALIZATION AND CALIBRATION
%==========================================================================
ncoils = length(coils2recon);
nslices = length(slices2recon);

for s = 1:nslices,
    slice = slices2recon(s);
    [filename,maxveloc,optr,nphases,nVE,nread,nintl,spiralid,...
        spatfov,spatres,pixels,kxkytraj,kxkyweights] = readdataparams(rawpath,slice);
    % %this loads the usercv variables from the rawdata file
    % filename = sprintf('%sslice%d.7',rawpath,slice);
    % [rawdata,usercv,hdr] = rawloadHD_jfn(filename,[],[],[], 1, [], [], []);
    % maxveloc = usercv(1); % maximum velocity value (1/2 of velocity field-of-view)
    % optr = usercv(2); % TR duration (in microseconds)
    % nphases = usercv(4); %number of cardiac phases (i.e., temporal frames)
    % nVE = usercv(7); %number of velocity encoding steps
    % nread = usercv(8); %number of readout samples
    % nintl = usercv(14); %number of spiral interleaves
    % spiralid = usercv(10); %this number identifies which spiral trajectory was used
    % heartrate = usercv(5); %heart rate of the subject during the scan
    % RRpct = usercv(6); %percentage of cardiac cycle that was covered
    % flipangle = usercv(15); %flip angle (in degrees)
    % % vesperbeat = usercv(3);
    % % vres = usercv(9);
    % % vegrad = usercv(11);
    % % realtime = usercv(12);
    % % nocine = usercv(13);
    % % ia_rf1 = usercv(16);
    % % densityreductionfactor = usercv(17);
    % % oblique = usercv(18);
    % % intlsperbeat = usercv(19);
    %
    % % load spatial parameters
    % switch(spiralid)
    %     case 19,
    %         kfile='recon16cm14mm8intl4vd';
    %         spatfov = 160;
    %         spatres = 1.4;
    %     otherwise,
    %         error(['unexpected spiralid value: ',num2str(spiralid)]);
    % end;
    % [kxkytraj kxkyweights] = kkread(kfile,nintl,nread);
    % % figure,plot(kxkytraj);axis equal; %plots kspace trajectory
    % kxkytraj = kxkytraj / (2*max(abs(kxkytraj(:)))); %normalize to [-0.5 , 0.5]
    % kxkyweights=kxkyweights/max(kxkyweights(:));%/sqrt(2); %normalizes the density values to [0 , sqrt(2)/2]
    % pixels = ceil(spatfov/spatres); %image size
    
    
    %==========================================================================
    % STARTING UP THE NUFFT ALGORITHM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('nufft_init...')
    
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
    xykvt = zeros(pixels,pixels,nVE,nphases); % dimensions are: x-y-kv-time
    res = zeros(pixels,pixels,nVE,nphases);
    for p = 1:nphases, %loops through the cardiac cycle
        xykvc = zeros(pixels,pixels,nVE,ncoils); % kv data in a each voxel for each coil
        rawdata = rawloadHD_jfn(filename,[],[],[], 1,1:ncoils,p,1:nintl); %reads one cardiac phase from each coil
        kxkyckv = permute(rawdata,[1 5 3 2 4]);
        for v = 1:nVE, %loops through the kv levels
            xyc = spirPI(kxkyckv(:,:,:,v),kxkyweights,kxkytraj,N,kSize,accel);
            xykvc(:,:,v,:) = permute(xyc,[1 2 4 3]); % = m.'; %flips the image and stores it
            disp(sprintf('kv=%2d/%d p=%3d/%d',v,nVE,p,nphases)) %display current cardiac phase
        end;
        
        xykvt(:,:,:,p) = combine4channels(xykvc(:,:,:,1),xykvc(:,:,:,2),xykvc(:,:,:,3),xykvc(:,:,:,4)); %combines the data from the 4 coils
        
    end;
    
    xyvt = fftshift(ifft(fftshift(xykvt,3),[],3),3); %the time velocity distsribution is calculated by inverse fourier transform along kv
    xy = xyvt(:,:,nVE/2+1,1); % final data, contains spatial image for v = 0 and t = 0
    
    %saves the data
    if(savedata2disk)
        disp('saving reconstructed data to disk...')
        cd datapaper\spirit
        save(sprintf('slice_%d_%d.mat',slice,accel), 'xyvt', 'xy', 'maxveloc', 'optr');
        disp('Done!')
        cd ..
    end
    
    figure,
    disp(' ')
    disp('Plotting data. Close the figure or press Ctrl-C to stop.')
    
    while(1),
        %plots the magnitude image (central kv, first cardiac phase)
        subplot(211),
        imshow(abs(xykvt(:,:,nVE/2+1,1)),[ ])
        set(gca,'YDir','normal')
        title(sprintf('slice %d',slice))
        
        %user will click on a pixel of the image
        xlabel('CLICK ON A BLOOD VESSEL / CLOSE FIGURE TO STOP')
        [Y,X] = ginput(1)
        
        %takes the velocity distribution from the prescribed pixel
        vt = xyvt(round(X),round(Y),:,:);
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
end