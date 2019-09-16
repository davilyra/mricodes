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
% Changes by Davi Marco Lyra-Leite <davi@ieee.org> and Joao Luiz Carvalho
% October 31, 2011

% More recent modification made by Davi M. Lyra-Leite, March 26 2012
% Inversion algorithm implemented (changing the position of FVE
% reconstruction steps - firstly it is implemented the iDFT in the kv axis,
% after the parallel imaging step)

%% starting the code by cleaning memory and closing opened figures
clear,clc,close all

% it will reconstruct data with undersampling factors of 1, 2 and 4, and
% will save the xyvt profile and xy image in a *.mat file
% if you want to work with only one slice coment the first loop, only one
% undersampling factor, coment the second and write its value
% slice = 4; % prescribe the slice number here (1 to 5)
% undersamplingfactor = 4;

usv = [1 2 4];
for slice = 1:5 % you can use slices from 1 to 5
    for us = 1:3
        undersamplingfactor = usv(us);

%==========================================================================
        rawpath = './rawdata/';
        nslices = 5;
        ncoils = 4;
        
%==========================================================================
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

%==========================================================================
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

        kxkytraj = kxkytraj / (2*max(abs(kxkytraj(:)))); %normalize to [-0.5 , 0.5]
        kxkyweights = kxkyweights/max(kxkyweights(:));%/sqrt(2); %normalizes the density values to [0 , sqrt(2)/2]


%==========================================================================
% starting up the SPIRiT image domain algorithm
%==========================================================================
        disp('spirit_init...')
        % Initial Parameters
        pixels = ceil(spatfov/spatres);      % image size
        kSize = [5,5];                       % SPIRiT Kernel size
        CalibSize = [30,30];                 % size of the calibration region
        N = [pixels,pixels];                 % size of the target image
        nIterCG = 50;                        % number of reconstruction iterations
        CalibTyk = 0.0125;                   % Tykhonov regularization for calibration 0.01-0.05 recommended
        lambda = 1;                          % Ratio between data and calibration consistency. 1 recommended when density compensation is used.
        accel = undersamplingfactor;         % Acceleration factor


%==========================================================================
        %reconstructs the time-velocity distributions from the selected slice
        disp(' '),disp(sprintf('reconstructing data from slice %d, undersampling factor %d...',slice,undersamplingfactor))
        filename = sprintf('%sslice%d.7',rawpath,slice); %filename of the current slice
        xyvt = zeros(pixels,pixels,nVE,nphases); % dimensions are: x-y-kv-time
        
        % loops through the cardiac cycle
        for p = 1:nphases
            xyvc = zeros(pixels,pixels,nVE,ncoils); % kv data in a each voxel for each coil
            xyvc1 = zeros(pixels,pixels,nVE,ncoils/2); % kv data in a each voxel for each coil
            rawdata = rawloadHD_jfn(filename,[],[],[], 1,2:3,p,1:nintl); % reads one cardiac phase from each coil
            rawdata2 = fftshift(ifft(fftshift(rawdata,2),[],2),2); % the time velocity distsribution is calculated by inverse fourier transform along kv
            kxkycv = permute(rawdata2,[1 5 3 2 4]);
            
            % loops through the v levels
            for v = 1:nVE
                xyc = spirPI(kxkycv(:,:,:,v),kxkyweights,kxkytraj,N,kSize,accel,CalibTyk,lambda);
                xyvc1(:,:,v,:) = permute(xyc,[1 2 4 3]);
                disp(sprintf('v =%3d/%d, p =%3d/%d',v,nVE,p,nphases)) % display current cardiac phase
                xyvc(:,:,v,1) = xyvc1(:,:,v,1);
                xyvc(:,:,v,2) = xyvc1(:,:,v,2);
            end
            
            xyvt(:,:,:,p) = combine4channels(xyvc(:,:,:,1),xyvc(:,:,:,2),xyvc(:,:,:,3),xyvc(:,:,:,4)); % combines the data from the 4 coils
        end
        
        xy = xyvt(:,:,nVE/2+1,1); % final data, contains spatial image for v = 0 and t = 0
        

%==========================================================================
        %saves the data
        disp('saving reconstructed data to disk...')
        cd datapaper\spirit_coil_selection
        save(sprintf('slice_%d_spirit_coils_23_%d.mat',slice,accel), 'xyvt', 'xy', 'maxveloc', 'optr');
        disp('Done!')
        cd ..
        cd ..

        
%==========================================================================
% plotting data
% if you want to select a voxel in order to see its velcotiy distribution,
% uncoment the while loop and coment the initial for loop which process the
% reconstruction over different undersampling factors and slices - so that
% it will be necessary to state which slice you want to recon and which
% undersampling factor you want to use

        figure,
        disp(' ')
        % disp('Plotting data. Close the figure or press Ctrl-C to stop.')

        % while(1),
        %    %plots the magnitude image (central kv, first cardiac phase)
        %     subplot(211),
            imshow(abs(xyvt(:,:,nVE/2+1,1)),[ ])
            set(gca,'YDir','normal')
            title(sprintf('slice %d,undersamplingfactor %d',slice,undersamplingfactor))
            pause(.1)

        %     %user will click on a pixel of the image
        %     xlabel('CLICK ON A BLOOD VESSEL / CLOSE FIGURE TO STOP')
        %     [Y,X] = ginput(1)
        % 
        %     %takes the velocity distribution from the prescribed pixel
        %     vt = xyvt(round(X),round(Y),:,:);
        %     vt = permute(vt,[3 4 1 2]);
        % 
        %     %plots the time-velocity distribution
        %     taxis = (0:(nphases-1))*optr/1000;
        %     vaxis = (((0:(nVE-1))/nVE)-0.5)*maxveloc*2;
        %     subplot(212),imagesc(taxis,vaxis,abs(vt))
        %     colormap(gray)
        %     set(gca,'YDir','normal')
        %     xlabel('time (ms)')
        %     ylabel('velocity (cm/s)')
        % 
        % end
    end
end