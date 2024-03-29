% RECONSTRUCTION OF THE TIME-VELOCITY DISTRIBUTION FROM
% MULTIPLE VOXELS IN A SLICE
%
% MATLAB WILL USE APPROXIMATELY 900MB OF RAM.
%
% Written by Joao L. A. Carvalho <joaoluiz@gmail.com>
% Department of Electrical Engineering
% University of Brasilia, Brazil
%
% July 14, 2008
%
% Modifications by Davi M. Lyra-Leite <davi@ieee.org>
% March 24, 2012

%% starting the code by cleaning memory and closing opened figures
clear,clc,close all


% it will reconstruct data with undersampling factors of 1, 2 and 4, and
% will save the xyvt profile and xy image in a *.mat file
% if you want to work with only one slice coment the first loop, only one
% undersampling factor, coment the second and write its value, eg.
% slice = 4; % prescribe the slice number here (1 to 5)
% undersamplingfactor = 4;

usv = [1 2 4];
for slice = 1:5
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

        kxkytraj = kxkytraj(:,1:undersamplingfactor:nintl);
        kxkyweights = kxkyweights(:,1:undersamplingfactor:nintl);

        kxkytraj = kxkytraj / (2*max(abs(kxkytraj(:)))); %normalize to [-0.5 , 0.5]
        kxkyweights=kxkyweights/max(kxkyweights(:));%/sqrt(2); %normalizes the density values to [0 , sqrt(2)/2]


%==========================================================================
% starting up the nuFFT algorithm
%==========================================================================
        disp('nufft_init...')
        pixels = ceil(spatfov/spatres); %image size
        Nd = [pixels,pixels];
        Jd = [6,6]; %kernel size (?)
        overgridfactor = 2;
        om = zeros(length(kxkytraj(:)),2);
        om(:,1) = 2*pi*real(kxkytraj(:)); %kspace trajectory
        om(:,2) = 2*pi*imag(kxkytraj(:));
        st = nufft_init(om, Nd, Jd, overgridfactor*Nd, Nd/2); %initialization
        
        
%==========================================================================
        % reconstructs the time-velocity distributions from the selected slice
        disp(' '),disp(sprintf('reconstructing data from slice %d, undersampling factor %d...',slice,undersamplingfactor))
        filename = sprintf('%sslice%d.7',rawpath,slice); % filename of the current slice
        xyvt = zeros(pixels,pixels,nVE,nphases); % dimensions are: x-y-kv-time
        
        % loops through the cardiac cycle
        for p = 1 : nphases
            xyvc = zeros(pixels,pixels,nVE,ncoils); % kv data in a each voxel for each coil
            
            % loops through each coil element
            for coil = 1:ncoils
                rawdata = rawloadHD_jfn(filename,[],[],[],1,coil,p,1:undersamplingfactor:nintl); % reads one cardiac phase from each coil
                rawdata2 = fftshift(ifft(fftshift(rawdata,2),[],2),2); % the time velocity distsribution is calculated by inverse fourier transform along kv
                
                % loops through the v levels
                for v = 1 : nVE
                    kxkydata = permute(rawdata2(:,v,1,1,:),[1 5 2 4 3]); % rearranges data from each kv level
                    weighteddata = kxkydata(:).*kxkyweights(:); % weights the kspace (kx-ky) data based on sampling density
                    m = nufft_adj(weighteddata,st)/pixels; % non-cartesian inverse Fourier transform along kx-ky
                    xyvc(:,:,v,coil) = m;
                end
                
            end
            
            xyvt(:,:,:,p) = (undersamplingfactor)*combine4channels(xyvc(:,:,:,1),xyvc(:,:,:,2),xyvc(:,:,:,3),xyvc(:,:,:,4)); % combines the data from the 4 coils
            disp(sprintf('p =%3d/%d',p,nphases)) % display current cardiac phase
        end;

        xy = xyvt(:,:,nVE/2+1,1); % final data, contains spatial image for v = 0 and t = 0

        
%==========================================================================
        %saves the data
        disp('saving reconstructed data to disk...')
        cd datapaper\sosnufft_3
        save(sprintf('slice_%d_sosnufft_%d.mat',slice,undersamplingfactor), 'xyvt','xy','maxveloc', 'optr');
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
%         disp('Plotting data. Close the figure or press Ctrl-C to stop.')
% 
%         while(1),

%            %plots the magnitude image (central kv, first cardiac phase)
%             subplot(211),
            imshow(abs(xyvt(:,:,nVE/2+1,1)),[ ])
            set(gca,'YDir','normal')
            title(sprintf('slice %d,undersamplingfactor %d',slice,undersamplingfactor))
            pause(.1)

%             %user will click on a pixel of the image
%             xlabel('CLICK ON A BLOOD VESSEL / CLOSE FIGURE TO STOP')
%             [Y,X] = ginput(1);X,Y,
%         
%             %takes the velocity distribution from the prescribed pixel
%             vt = xyvt(round(X),round(Y),:,:);
%             vt = permute(vt,[3 4 1 2]);
%         
%             %plots the time-velocity distribution
%             taxis = (0:(nphases-1))*optr/1000;
%             vaxis = (((0:(nVE-1))/nVE)-0.5)*maxveloc*2;
%             subplot(212),imagesc(taxis,vaxis,abs(vt))
%             colormap(gray)
%             set(gca,'YDir','normal')
%             xlabel('time (ms)')
%             ylabel('velocity (cm/s)')

%         end
    end
end
