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

% New modifications by Joao L. A. Carvalho <joaoluiz@gmail.com>, July 17
% 2012 with new view-ordering algorithm, implementation of a single
% calibration step based on time-mean of undersampled data in order to
% generate the fully sampled pattern in the center


%% starting the code by cleaning memory and closing opened figures
clear,clc,close all

% it will reconstruct data with undersampling factors of 1, 2 and 4, and
% will save the xyvt profile and xy image in a *.mat file
% if you want to work with only one slice coment the first loop, only one
% undersampling factor, coment the second and write its value
% slice = 4; % prescribe the slice number here (1 to 5)
% undersamplingfactor = 4;

%==========================================================================
rawpath = './rawdata/';
nslices = 5;
ncoils = 4;

%==========================================================================

usv = [1 2 4];
% usv = 4;
for slice = 4, %1:5 % you can use slices from 1 to 5
    for us = 1:length(usv),
        undersamplingfactor = usv(us);

        %this loads the usercv variables from the rawdata file
        filename = sprintf('%sslice%d.7',rawpath,slice); %filename of the current slice
        [~,usercv,hdr] = rawloadHD_jfn(filename,[],[],[], 1, [], [], []);
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

        % nread = 1012;
        kxkytraj = kxkytraj(1:nread,:);
        kxkyweights = kxkyweights(1:nread,:);
        kxkytraj = kxkytraj*256/spatfov; %corrects k-space trajectory
        spatres = 1/(2*abs(kxkytraj(end,1))); %calculates spatial resolution from kspace trajectory
        %figure,plot(kxkytraj);axis equal; %plots kspace trajectory
        kxkytraj = kxkytraj*spatres; %normalize to [-0.5 , 0.5]
        kxkyweights = kxkyweights/max(kxkyweights(:));%/sqrt(2); %normalizes the density values to [0 , sqrt(2)/2]


%==========================================================================
% starting up the SPIRiT image domain algorithm
%==========================================================================
        disp('spirit_init...')
        % Initial Parameters
        pixels = ceil(spatfov/spatres);     % image size
        kSize = [5,5];                      % SPIRiT Kernel size
        CalibSize = [30,30];                % size of the calibration region
        imagesize = [pixels,pixels];        % size of the target image
        CalibTyk = 0.0125;                  % Tykhonov regularization for calibration 0.01-0.05 recommended
        lambda = 1;
        nIterCG = 50;                       % number of reconstruction iterations


%==========================================================================
        %reconstructs the time-velocity distributions from the selected slice
        fprintf('reconstructing data from slice %d, undersampling factor %d...',slice,undersamplingfactor)
        disp(' ')


%==========================================================================
% SPIRiT Calibration

        disp(' '),fprintf('starting calibration...'), disp(' ')

        phasesCalib = 1:floor(nphases/undersamplingfactor)*undersamplingfactor; %(nphases-undersamplingfactor+1):nphases;
        nphasesCalib = length(phasesCalib);
        kvCalib = nVE/2 + 1;
        rawdata = rawloadHD_jfn(filename,[],[],[], 1,1:ncoils,phasesCalib,1:nintl); % reads one cardiac phase from each coil
        %size(rawdata) = [readout kv coil phases intls]

        % undersampling pattern
        rkvcpi = nan(nread,1,ncoils,nphasesCalib/undersamplingfactor,nintl);
        for k = 1%1:undersamplingfactor
            idx = (k:undersamplingfactor:nintl);
            phases2 = k:undersamplingfactor:nphasesCalib;
            rkvcpi(:,1,:,:,idx) = rawdata(1:nread,kvCalib,:,phases2,idx);
        end
        kxkyc = permute(mean(rkvcpi,4),[1 5 3 2 4]);
    
        % Gridding and density compensation reconstruction
        GFFT1 = NUFFT(kxkytraj,kxkyweights, [0,0] , imagesize);
        im = GFFT1'*(kxkyc.*repmat(sqrt(kxkyweights),[1,1,ncoils]));
        %imagesc(abs(im(:,:,1))),colorbar,error
        kData = fft2c(im);
        kernel = zeros([kSize,ncoils,ncoils]);
        kCalib = crop(kData,[CalibSize,ncoils]); % crop center k-space for calibration

        % Calibration matrix (more information, refeers to Phill Beatty's
        % work)
        [AtA,] = corrMatrix(kCalib,kSize);
        for n=1:ncoils
            kernel(:,:,:,n) = calibrate(AtA,kSize,ncoils,n,CalibTyk);
        end
        GOP = SPIRiT(kernel, 'image',imagesize); % init the SPIRiT Operator


%==========================================================================
% SPIRiT Reconstruction

        disp(' '),fprintf('starting SPIRiT reconstruction...'),disp(' ')

        % loops through the cardiac cycle
        xyvt = zeros(pixels,pixels,nVE,nphases); % dimensions are: x-y-kv-time
        for p = 1:nphases
            fprintf('p =%3d/%d',p,nphases) % display current cardiac phase
            disp(' ')
            
            intls2recon = (mod(p - 1,undersamplingfactor) + 1):undersamplingfactor:nintl;
            w_u = kxkyweights(:,intls2recon);
            k_u = kxkytraj(:,intls2recon);
            GFFT_u = NUFFT(k_u,w_u, [0,0], imagesize);
            wu2 = repmat(sqrt(w_u),[1,1,ncoils]);
            
            rawdata = rawloadHD_jfn(filename,[],[],[], 1,1:ncoils,p,intls2recon); % reads one cardiac phase from each coil
            
            % the time velocity distsribution is calculated by inverse fourier transform along kv
            kxkycv = undersamplingfactor*permute(fftshift(ifft(fftshift(rawdata(1:nread,:,:,:,:),2),[],2),2),[1 5 3 2 4]); 
%             kxkycv(:,:,3,:) = 0;
%             kxkycv(:,:,4,:) = 0;
            
            % loops through the velocity values
            xyvc = zeros(pixels,pixels,nVE,ncoils); % kv data in a each voxel for each coil
            %tic
            for v = 1:nVE,
%             parfor v = 1:nVE, %PERGUNTAR PARA A ROSANA PORQUE NAO ESTA ACELERANDO
                xyvc(:,:,v,:) = permute(cgNUSPIRiT(kxkycv(:,:,:,v).*wu2,GFFT_u'*(kxkycv(:,:,:,v).*wu2),GFFT_u,GOP,nIterCG,lambda),[1 2 4 3]);
                %disp(sprintf('v =%3d/%d, p =%3d/%d',v,nVE,p,nphases)) % display current cardiac phase

                % %parfor variables
                % wu2,GFFT_u,GOP,nIterCG,lambda: broadcast variables
                % v: loop variable
                % kxkycv: input sliced variable
                % xyvc: output sliced variable
            end
            %toc,error

            xyvt(:,:,:,p) = sqrt( xyvc(:,:,:,1).^2 + xyvc(:,:,:,2).^2 + xyvc(:,:,:,3).^2 + xyvc(:,:,:,4).^2); % combines the data from the 4 coils
%             xyvt(:,:,:,p) = combine4channels( xyvc(:,:,:,1),xyvc(:,:,:,2),xyvc(:,:,:,3),xyvc(:,:,:,4)); % combines the data from the 4 coils
        end

        xy = xyvt(:,:,nVE/2+1,1); % final data, contains spatial image for v = 0 and t = 0
        

%==========================================================================
%saves the data

        disp('saving reconstructed data to disk...')
        save(sprintf('datapaper/spirit_JL1/recon_slice_%d_%d.mat',slice,undersamplingfactor), 'xyvt', 'xy', 'maxveloc', 'optr');
%         save(sprintf('datapaper/spirit_JL1/recon_slice_%d_spirit_coils_12_%d.mat',slice,undersamplingfactor), 'xyvt', 'xy', 'maxveloc', 'optr');
        disp('Done!')

        
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