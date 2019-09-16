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
%
% Current version by Davi M. Lyra-Leite: Apr 22, 2013


%% starting the code by cleaning memory and closing opened figures
clear,clc,close all

%% setting code parameters
savedata2disk = 0; % [0] no, [1] yes
rawpath = './rawdata/';
slices2recon = [1]; %examples: [4] or [1 2 3 4 5]
coils2recon = [1 2 3 4]; %[1 2 3 4]


%==========================================================================
% INITIALIZATION AND CALIBRATION
%==========================================================================
ncoils = length(coils2recon);
nslices = length(slices2recon);

for s = 1:nslices,
    slice = slices2recon(s);
    [filename,maxveloc,optr,nphases,nVE,nread,nintl,spiralid,...
        spatfov,spatres,pixels,kxkytraj,kxkyweights] = readdataparams(rawpath,slice);
        
    
%==========================================================================
% STARTING UP THE NUFFT ALGORITHM
%==========================================================================
    disp('nufft_init...')
    
    Nd = [pixels,pixels];
    Jd = [6,6]; %kernel size (?)
    overgridfactor = 2;
    om(:,1) = 2*pi*real(kxkytraj(:)); %kspace trajectory
    om(:,2) = 2*pi*imag(kxkytraj(:));
    st = nufft_init(om, Nd, Jd, overgridfactor*Nd, Nd/2); %initialization
    
    
%==========================================================================
% STARTING UP THE SPIRiT ALGORITHM
%==========================================================================
    % Initial Parameters
    kSize = [5,5];                    % SPIRiT Kernel size
    CalibSize = [30,30];              % size of the calibration region
    N = Nd;                           % size of the target image
    nIterCG = 50;                     % number of reconstruction iterations
    CalibTyk = 0.02;                  % Tykhonov regularization for calibration 0.01-0.05 recommended
    lambda = 1;                       % Ratio between data and calibration consistency. 1 recommended when density compensation is used.
    accel = 4;                        % Acceleration factor
    

%==========================================================================
% RECONSTRUCTION
%==========================================================================
    
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