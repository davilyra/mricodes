% RECONSTRUCTION OF THE TIME-VELOCITY DISTRIBUTION FROM
% MULTIPLE VOXELS IN A SLICE, WHICH ARE MANUALLY
% PRESCRIBED BY THE USER. MATLAB WILL USE APPROXIMATELY 1.1 GB OF RAM.
%
% This code allows the use of 2D SPIRiT or nuFFT+SoS as reconstruction
% methods. It also enables the user to select between the original
% reconstruction of spiral FVE data (as presentend by Carvalho and Nayak in
% MRM 57:639, 2007) or the proposed framework presented by Lyra-Leite et
% al. in EMBC 34:416.
%
% Original reconstruction algorithm by  Joao Carvalho <joaoluiz@gmail.com>
% Department of Electrical Engineering
% University of Brasilia, Brazil
% July 14, 2008
%
% Current algorithm's version by Davi M. Lyra-Leite
% Department of Electrical Engineering
% University of Southern California, Los Angeles
% United States
% Apr 22, 2013


%% starting the code by cleaning memory and closing opened figures
clear all
clc
close all

%% setting code parameters
savedata2disk = 1; % [0] no, [1] yes
rawpath = './rawdata/';

displayvfspace = 1; % this is will increase reconstruction time, because it requires combinining the phase images
dbrange = 40; % dynamic range (in db) for displaying the v-f space (recommended value: 40). Use dbrange = 0 for linear scale.

slices2recon = [4]; % examples: [4] or [1 2 3 4 5]
coils2recon = [1 2 3 4]; % [1 2 3 4]
accelfactors = [2]; % [1 2 4];
useinversion = 0; % [0] no, [1] yes - this follows what is presented in Lyra-Leite et al. EMBC 34:416, 2012.
usespirit = 1; % set to 1 to use spirit and 0 to use SoS+nuFFT
plotdata = 1; % set 1 to plot the reconstructed image, 0 otherwise

%==========================================================================
% INITIALIZATION AND CALIBRATION
%==========================================================================
ncoils = length(coils2recon);
nslices = length(slices2recon);


if usespirit == 1
%==========================================================================
% STARTING UP THE SPIRiT ALGORITHM
%==========================================================================
    % Initial Parameters
    kSize = [5,5];                    % SPIRiT Kernel size
    CalibSize = [30,30];              % size of the calibration region
    nIterCG = 50;                     % number of reconstruction iterations
    CalibTyk = 0.02;                  % Tykhonov regularization for calibration 0.01-0.05 recommended
    lambda = 1;                       % Ratio between data and calibration consistency. 1 recommended when density compensation is used.
end


for s = 1:nslices,
    slice = slices2recon(s);
    
    for us = 1:length(accelfactors),
        
        undersamplingfactor = accelfactors(us);
        [filename,maxveloc,optr,nphases,nVE,nread,nintl,spiralid,...
            spatfov,spatres,pixels,kxkytraj,kxkyweights] = readdataparams(rawpath,slice);

        if ~(usespirit)
            fprintf('starting NUFFT+SOS reconstruction...')
        else
            fprintf('starting NUFFT+SPIRiT+SOS reconstruction...')
        end
        [GFFT,wroot,intls2read] = initrecon(ncoils,undersamplingfactor,kxkytraj,kxkyweights,pixels,pixels);
        [kmatrix,~,~] = calculatekmatrix(nVE,nphases,0,undersamplingfactor);
        
        
%==========================================================================
% RECONSTRUCTION
%==========================================================================
        
        %reconstructs the time-velocity distributions from the selected slice
        disp(' '),fprintf('reconstructing data from slice %d...',slice)
        disp(' '),fprintf('undersamplingfactor of %d...',undersamplingfactor)
        
        xykvt = zeros(pixels,pixels,nVE,nphases); % dimensions are: x-y-kv-time
        xyvt = zeros(pixels,pixels,nVE,nphases); % dimensions are: x-y-v-time
        
        for p = 1:nphases, %loops through the cardiac cycle
%         for p = 1:1, %loops through the cardiac cycle
            rawdata = rawloadHD_jfn(filename,[],[],[], 1,coils2recon,p,1:nintl); %reads one cardiac phase from each coil
            rawdata2recon = zeros(nread,nVE,ncoils,1,nintl);
            kweights2recon = zeros(nread,nintl);
            ktraj2recon = zeros(nread,nintl);
            
            for v = 1:nVE,
                k = kmatrix(v,p);
                if ~isnan(k),
                    rawdata2recon(:,v,:,:,intls2read(:,k)) = rawdata(1:nread,v,:,:,intls2read(:,k));
                    kweights2recon(:,intls2read(:,k)) = kxkyweights(1:nread,intls2read(:,k));
                    ktraj2recon(:,intls2read(:,k)) = kxkytraj(1:nread,intls2read(:,k));
                end
            end
            
            xykvc = zeros(pixels,pixels,nVE,ncoils); % kv data in a each voxel for each coil
            xyvc = zeros(pixels,pixels,nVE,ncoils); % kv data in a each voxel for each coil
            
            if useinversion == 0
                clear xyvc
                kxkyckv = permute(rawdata2recon,[1 5 3 2 4]);
                
                if ~(usespirit)
                    disp(' ')
                    fprintf('p=%3d/%d',p,nphases) %display current cardiac phase
                    
                    for v = 1:nVE, %loops through the kv levels
                        xykvc(:,:,v,:) = undersamplingfactor*permute(GFFT'*(wroot.*kxkyckv(:,:,:,v)),[1 2 4 3]);
                    end
                    
                    xykvt(:,:,:,p) = combine_coils(xykvc,coils2recon);
                    xyvt = fftshift(ifft(fftshift(xykvt,3),[],3),3); %the time velocity distsribution is calculated by inverse fourier transform along kv
                    
                else
                    kxkyckv_fullysampled = permute(rawdata,[1 5 3 2 4]);
                    
                    for v = 1:nVE, %loops through the kv levels
                        xyc = spirPI3(kxkyckv_fullysampled(:,:,:,v),kxkyckv(:,:,:,v),kxkyweights,kweights2recon,kxkytraj,ktraj2recon,[pixels,pixels],kSize,undersamplingfactor);
                        xykvc(:,:,v,:) = permute(xyc,[1 2 4 3]);
                        disp(' ')
                        fprintf('kv=%3d/%d, p=%3d/%d',v,nVE,p,nphases) %display current cardiac phase
                    end
                    
                    xykvt(:,:,:,p) = combine_coils(xykvc,coils2recon);
                    xyvt = fftshift(ifft(fftshift(xykvt,3),[],3),3); %the time velocity distsribution is calculated by inverse fourier transform along kv
                end
                
            else
                clear xykvc
                rawdata2 = fftshift(ifft(fftshift(rawdata2recon,2),[],2),2); % the time velocity distsribution is calculated by inverse fourier transform along kv
                kxkycv = permute(rawdata2,[1 5 3 2 4]);
                
                if ~(usespirit)
                    disp(' ')
                    fprintf('p =%3d/%d',p,nphases) % display current cardiac phase
                    for v = 1:nVE,
                        xyvc(:,:,v,:) = undersamplingfactor*permute(GFFT'*(wroot.*kxkycv(:,:,:,v)),[1 2 4 3]);
                    end
                    xyvt(:,:,:,p) = combine_coils(xyvc,coils2recon);
                    
                else
                    for v = 1:nVE % loops through the v levels
                        xyc = spirPI(kxkycv(:,:,:,v),kxkyweights,kxkytraj,[pixels,pixels],kSize,undersamplingfactor,CalibTyk,lambda);
                        xyvc(:,:,v,:) = permute(xyc,[1 2 4 3]);
                        disp(' ')
                        fprintf('v =%3d/%d, p =%3d/%d',v,nVE,p,nphases) % display current cardiac phase
                    end
                    xyvt(:,:,:,p) = combine_coils(xyvc,coils2recon);
                end
                
            end
            
        end
        xy = xyvt(:,:,nVE/2+1,1); % final data, contains spatial image for v = 0 and t = 0
        
        if (savedata2disk) && ncoils == 4
            disp(' ')
            disp('saving reconstructed data to disk...')
            if ~(usespirit)
                if ~(useinversion)
                    save(sprintf('recondata/spirit_599/recon_sos_no_inversion_slice_%d_us_%d.mat',slice,undersamplingfactor), 'xyvt', 'xy', 'maxveloc', 'optr');
                else
                    save(sprintf('recondata/spirit_599/recon_sos_inversion_slice_%d_us_%d.mat',slice,undersamplingfactor), 'xyvt', 'xy', 'maxveloc', 'optr');
                end
            else
                if ~(useinversion)
                    save(sprintf('recondata/spirit_599/recon_l1spirit_no_inversion_slice_%d_us_%d.mat',slice,undersamplingfactor), 'xyvt', 'xy', 'maxveloc', 'optr');
                else
                    save(sprintf('recondata/spirit_599/recon_l1spirit_inversion_slice_%d_us_%d.mat',slice,undersamplingfactor), 'xyvt', 'xy', 'maxveloc', 'optr');
                end
            end
            disp('Done!')
            disp(' ')
        end
        
        if plotdata == 1
            %plots the magnitude image (v = 0, first cardiac phase)
            figure,
            subplot(211),
            imshow(abs(xy),[ ])
            set(gca,'YDir','normal')
            title(sprintf('slice %d,undersamplingfactor %d',slice,undersamplingfactor))
            
            disp(' ')
            disp('Plotting data. Close the figure or press Ctrl-C to stop.')
            while(1),
                %user will click on a pixel of the image
                [Y,X] = ginput(1);
                vt = permute(xyvt(round(X),round(Y),:,:),[3 4 1 2]); %takes the velocity distribution from the prescribed pixel
                displayvtvf(vt,optr,maxveloc,dbrange,X,Y,displayvfspace,slice)
            end
        end
    end
end