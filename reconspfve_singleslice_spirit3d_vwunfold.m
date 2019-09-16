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

% New modifications by Joao L. A. Carvalho <joaoluiz@gmail.com>, July 31
% 2012 with new view-ordering algorithm, implementation of a single
% calibration step based on time-mean of undersampled data in order to
% generate the fully sampled pattern in the center

% Reconstruction uses the 3D Spiral SPIRiT, developed by Taehoon Shin, and
% modified by Davi M. Lyra-Leite. First version: September 26th 2012.

% Current version: Apr 21st 2013.


%% starting the code by cleaning memory and closing opened figures
clear, clc, close all


%% setting code parameters
displayvfspace = 1; % this is will increase reconstruction time, because it requires combinining the phase images
dbrange = 40; % dynamic range (in db) for displaying the v-f space (recommended value: 40). Use dbrange = 0 for linear scale.

%==========================================================================
savedata2disk = 0; % [0] no, [1] yes
rawpath = './rawdata/';
slices2recon = [1]; %examples: [4] or [1 2 3 4 5]
voxel2recon = [ ]; %[X,Y], or [ ] to select from image in while loop (CCA: [74,71])
coils2recon = [1 2 3 4]; %[1 2 3 4]
accelfactors = [4]; %[1 2 4];

viewordering = 4;
% [0] same interleaves for all phases (aliasing will be at v = 0 and f =0)
% [1] alternate interleaves between phases (interleaf view-sharing) (aliasing will be at v = 0 and shifted in f)
% [2] alternate interleaves between kvs (aliasing will be at f = 0 and shifted in v)
% [3] alternate interleaves between phases and kvs (spatial aliasing will be shifted in v and f)
% [4] unfold-style view ordering (aliasing will be shifted in v and f)
% [11,12,13]: similar to 1, 2 and 3 but with worse results
% [21,22,23]: similar to 1, 2 and 3 but with random sampling


%==========================================================================
% INITIALIZATION AND CALIBRATION
%==========================================================================

% initialization;

%==========================================================================
nslices = length(slices2recon);
ncoils = length(coils2recon);
%==========================================================================


for s = 1:nslices,
    slice = slices2recon(s);
    
    for us = 1:length(accelfactors),
        undersamplingfactor = accelfactors(us);
        
        [filename,maxveloc,optr,nphases,nVE,nread,nintl,spiralid,...
            spatfov,spatres,pixels,kxkytraj,kxkyweights] = readdataparams(rawpath,slice);

        [kmatrix,iuf,vuf] = calculatekmatrix(nVE,nphases,viewordering,undersamplingfactor);

        fprintf('reconstructing data from slice %d, undersampling factor %d...',slice,undersamplingfactor)

        % SPIRiT Calibration single step
        imSize = [pixels pixels];
        [GOP, sampidx_tot, FT_SPR3d, f0] = spiritcalibration3d_davi(kmatrix,coils2recon,undersamplingfactor,iuf,filename,kxkytraj,kxkyweights,imSize);
        lambda = 2;
        nIterCG = 15;
        disp(' ')
        fprintf('starting NUFFT+SPIRiT+SOS reconstruction...')
        disp(' ')
        
        [GFFT,wroot,intls2read,w_u] = initrecon3(ncoils,iuf,kxkytraj,kxkyweights,pixels,pixels,nVE);
        wroot = permute(wroot,[1 2 4 3]);
        

%==========================================================================
% RECONSTRUCTION
%==========================================================================
        % loops through the cardiac cycle
        xyvt = zeros(pixels,pixels,nVE,nphases);
        for p = 1:nphases,
            fprintf('p =%3d/%d',p,nphases) % displays current cardiac phase number
            disp(' ')
            rawdatakv = rawloadHD_jfn(filename,[],[],[],1,coils2recon,p,1:nintl); % reads one cardiac phase: %read,kv,coils,phases,intl
            rawdata2recon = zeros(nread,nVE,ncoils,1,nintl);
            
            for v = 1:nVE,
                k = kmatrix(v,p);
                if ~isnan(k),
                    rawdata2recon(:,v,:,:,intls2read(:,k)) = rawdatakv(1:nread,v,:,:,intls2read(:,k));
                end
            end            
%             rawdatav = fftshift(ifft(fftshift(rawdata2recon,2),[],2),2); % the time velocity distsribution is calculated by inverse fourier transform along kv
            kxkyckv = permute(rawdata2recon,[1 5 2 3 4]); % size: nread, intl, nVE, coils

            GFFT_u = SPR3dFULL(GFFT, [ pixels,pixels, nVE], [nread,size(intls2read),nVE], ncoils, sampidx_tot, f0 ); % SPIRiT 3D Operator

            imund = GFFT_u'*(kxkyckv.*repmat(sqrt(w_u),[1,1,nVE,ncoils]))*undersamplingfactor; % undersampled image
            
%             xykvc = zeros(pixels,pixels,nVE,ncoils); % kv data in each voxel for each coil
            tic;
            [ xykvc FLAG,~,~,~,~]= cgNUSPIRiT3d(double(kxkyckv).*wroot,imund,FT_SPR3d,GOP,nIterCG,lambda);
            toc;
            

%==========================================================================
% SUM OF SQUARES (combines the data from the different coils)
%==========================================================================
            for c = 1:ncoils,
                xyvt(:,:,:,p) = xyvt(:,:,:,p) + abs(xykvc(:,:,:,c)).^2;
            end
            xyvt(:,:,:,p) = sqrt(xyvt(:,:,:,p));
        end
        xy = xyvt(:,:,floor(nVE/2+1),1); % spatial image for v = 0 and t = 0
        disp(' ')

%==========================================================================
% DISPLAY AND SAVE TO DISK
%==========================================================================
        if(savedata2disk)
            disp('saving reconstructed data to disk...')
            save(sprintf('datapaper/spirit_JL1/recon_spirit_slice_%d_us_%d_method_%d.mat',slice,undersamplingfactor,viewordering), 'xyvt', 'xy', 'maxveloc', 'optr');
            disp('Done!')
            disp(' ')
        end

        %plots the magnitude image (v = 0, first cardiac phase)
        figure(viewordering + 1),subplot(211),
        imshow(abs(xy),[ ])
        set(gca,'YDir','normal')
        title(sprintf('slice %d,undersamplingfactor %d, viewordering %d',slice,undersamplingfactor,viewordering))

        if isempty(voxel2recon)
            disp(' ')
            disp('Plotting data. Close the figure or press Ctrl-C to stop.')
            xlabel('CLICK ON A BLOOD VESSEL / CLOSE FIGURE TO STOP')
            while(1),
                %user will click on a pixel of the image
                [Y,X] = ginput(1);
                vt = permute(xyvt(round(X),round(Y),:,:),[3 4 1 2]); %takes the velocity distribution from the prescribed pixel
                displayvtvf(vt,optr,maxveloc,dbrange,X,Y,displayvfspace,slice)
            end
        else
            X = voxel2recon(1);
            Y = voxel2recon(2);
            vt = permute(xyvt(round(X),round(Y),:,:),[3 4 1 2]); %takes the velocity distribution from the prescribed pixel
            displayvtvf(vt,optr,maxveloc,dbrange,X,Y,displayvfspace,slice)
        end
    end
end