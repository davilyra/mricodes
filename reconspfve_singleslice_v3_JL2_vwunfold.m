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


%% starting the code by cleaning memory and closing opened figures
clear, clc, close all


%% setting code parameters
usespirit = 0; % set to 1 to use spirit, 0 otherwise
useparfor = 0; % set to 1 to use parfor in multi-core CPUs (implemented only for spirit reconstruction)

displayvfspace = 1; % this is will increase reconstruction time, because it requires combinining the phase images
reconpsf = 0; % set to 1 to reconstruct the point spread function instead of the actual data (works only with usespirit = 0)
dbrange = 40; % dynamic range (in db) for displaying the v-f space (recommended value: 40). Use dbrange = 0 for linear scale.

%==========================================================================
savedata2disk = 1;
rawpath = './rawdata/';
slices2recon = [4]; %examples: [4] or [1 2 3 4 5]
voxel2recon = [74,71]; %[X,Y], or [ ] to select from image in while loop (CCA: [74,71])
coils2recon = [1 2 3 4]; %[1 2 3 4]
accelfactors = [2 4]; %[1 2 4];

viewordering = 21;
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
            spatfov,spatres,pixels,kxkytraj,kxkyweights]=readdataparams(rawpath,slice);

        [kmatrix,iuf,vuf] = calculatekmatrix(nVE,nphases,viewordering,undersamplingfactor);
        
        %kmatrix,error
        
        disp(' '),disp(sprintf('reconstructing data from slice %d, undersampling factor %d...',slice,undersamplingfactor))

        if(usespirit), % SPIRIT CALIBRATION
            [GOP,lambda,nIterCG] = spiritcalibration2(kmatrix,coils2recon,undersamplingfactor,iuf,nread,nintl,filename,kxkytraj,kxkyweights,pixels,pixels);
            disp(' '),disp(sprintf('starting NUFFT+SPIRiT+SOS reconstruction...'))
        else
            disp(' '),disp(sprintf('starting NUFFT+SOS reconstruction...'))
        end;
        [GFFT,wroot,intls2read] = initrecon(ncoils,iuf,kxkytraj,kxkyweights,pixels,pixels);

        if(reconpsf)
            nphases = 4;% undersamplingfactor;
            nVE = 4; %undersamplingfactor;
            kmatrix = kmatrix(1:nVE,1:nphases);
%         else,
%             nphases = 4, nphases = 4, nphases = 4,
%             nphases = 4, nphases = 4, nphases = 4,
%             nphases = 4, nphases = 4, nphases = 1,
%             kmatrix = kmatrix(1:nVE,1:nphases),
        end;
        
        
%==========================================================================
% RECONSTRUCTION
%==========================================================================

        % loops through the cardiac cycle
        xyvt = zeros(pixels,pixels,nVE,nphases);
        for p = 1:nphases,
            disp(sprintf('p =%3d/%d',p,nphases)) % displays current cardiac phase number
            rawdatakv = rawloadHD_jfn(filename,[],[],[],1,coils2recon,p,1:nintl); % reads one cardiac phase: %read,kv,coils,phases,intl
            rawdata2recon = zeros(nread,nVE,ncoils,1,nintl);
            for v = 1:nVE,
                k = kmatrix(v,p);
                if ~isnan(k),
                    if reconpsf,
                        rawdata2recon(:,v,:,:,intls2read(:,k)) = ones([nread,1,ncoils,1,size(intls2read,1)]);
                    else
                        rawdata2recon(:,v,:,:,intls2read(:,k)) = rawdatakv(1:nread,v,:,:,intls2read(:,k));
                    end
                end
            end
            
            rawdatav = fftshift(ifft(fftshift(rawdata2recon,2),[],2),2); % the time velocity distsribution is calculated by inverse fourier transform along kv
            kxkycv = permute(rawdatav,[1 5 3 2 4]);

            % loops through the velocity values
            xyvc = zeros(pixels,pixels,nVE,ncoils); % kv data in a each voxel for each coil
            if usespirit,
                tic;
                if useparfor,
                    % list of parfor variables
                    % --------------------------
                    % v: loop variable
                    % kxkycv: sliced input variable
                    % xyvc: sliced output variable
                    % weightedraw,xyc: temporary variables
                    % wroot,GFFT,GOP,nIterCG,lambda: broadcast variables
                    parfor v = 1:nVE,
                        weightedraw = wroot.*kxkycv(:,:,:,v);
                        x0 = GFFT'*weightedraw;
                        xyc = cgNUSPIRiT(weightedraw,x0,GFFT,GOP,nIterCG,lambda); %NUFFT + SPIRIT
%                         xyc = cgNUSPIRiTJL(weightedraw,GFFT,GOP,nIterCG,lambda); %NUFFT + SPIRIT
                        xyvc(:,:,v,:) = permute(xyc,[1 2 4 3]);
                    end
                else
                    for v = 1:nVE,
                        weightedraw = wroot.*kxkycv(:,:,:,v);
                        x0 = GFFT'*weightedraw;
                        xyc = cgNUSPIRiT(wroot.*kxkycv(:,:,:,v),x0,GFFT,GOP,nIterCG,lambda); %NUFFT + SPIRIT
%                         xyc = cgNUSPIRiTJL(wroot.*kxkycv(:,:,:,v),GFFT,GOP,nIterCG,lambda); %NUFFT + SPIRIT
                        xyvc(:,:,v,:) = permute(xyc,[1 2 4 3]);
                        disp(sprintf('p =%3d/%d, v =%3d/%d',p,nphases,v,nVE)) % display current cardiac phase
                    end
                end
                toc;
            else %NUFFT only
                for v = 1:nVE,
                    % reconstructed image is scaled by "undersamplingfactor" to compensate
                    % for reduced ammount of signal due to data undersampling
                    xyvc(:,:,v,:) = undersamplingfactor*permute(GFFT'*(wroot.*kxkycv(:,:,:,v)),[1 2 4 3]);
                end
            end

            
%==========================================================================
            % SUM OF SQUARES (combines the data from the different coils)
            if reconpsf || displayvfspace,
                switch ncoils,
                    case 1,
                        xyvt(:,:,:,p) = xyvc(:,:,:,1);
                    case 2,
                        xyvt(:,:,:,p) = combine4channels(xyvc(:,:,:,1),xyvc(:,:,:,2),0,0);
                    case 3,
                        xyvt(:,:,:,p) = combine4channels(xyvc(:,:,:,1),xyvc(:,:,:,2),xyvc(:,:,:,3),0);
                    case 4,
                        xyvt(:,:,:,p) = combine4channels(xyvc(:,:,:,1),xyvc(:,:,:,2),xyvc(:,:,:,3),xyvc(:,:,:,4));
                end
            else
                for c=1:ncoils,
                    xyvt(:,:,:,p) = xyvt(:,:,:,p) + abs(xyvc(:,:,:,c)).^2;
                end
                xyvt(:,:,:,p) = sqrt(xyvt(:,:,:,p));
            end

        end

        xy = xyvt(:,:,floor(nVE/2+1),1); % spatial image for v = 0 and t = 0

        
%==========================================================================
% DISPLAY AND SAVE TO DISK
%==========================================================================

        if(reconpsf) %displays point-spread-function  
            showpsf(xyvt,maxveloc,optr,spatfov,dbrange);
        else
            if(savedata2disk)
                if ~(usespirit)
                    disp('saving reconstructed data to disk...')
                    cd datapaper/spirit_JL1
                    save(sprintf('recon_slice_sos_slice_%d_us_%d_method_%d.mat',slice,undersamplingfactor,viewordering), 'xyvt', 'xy', 'maxveloc', 'optr');
                    disp('Done!')
                    cd ..
                    cd ..
                else
                    disp('saving reconstructed data to disk...')
                    cd datapaper/spirit_JL1
                    save(sprintf('recon_slice_spirit_slice_%d_us_%d_method_%d.mat',slice,undersamplingfactor,viewordering), 'xyvt', 'xy', 'maxveloc', 'optr');
                    disp('Done!')
                    cd ..
                    cd ..
                end
            end

            %plots the magnitude image (v = 0, first cardiac phase)
            figure(viewordering),subplot(211),
            imshow(abs(xy),[ ])
            set(gca,'YDir','normal')
            title(sprintf('slice %d,undersamplingfactor %d',slice,undersamplingfactor))
            

            if isempty(voxel2recon)
                disp(' ')
                disp('Plotting data. Close the figure or press Ctrl-C to stop.')
                xlabel('CLICK ON A BLOOD VESSEL / CLOSE FIGURE TO STOP')
                while(1),
                    %user will click on a pixel of the image
                    [Y,X] = ginput(1);
                    vt = permute(xyvt(round(X),round(Y),:,:),[3 4 1 2]); %takes the velocity distribution from the prescribed pixel
                    displayvtvf(vt,optr,maxveloc,dbrange,X,Y,slice,displayvfspace)
                end
            else
                X = voxel2recon(1);
                Y = voxel2recon(2);
                vt = permute(xyvt(round(X),round(Y),:,:),[3 4 1 2]); %takes the velocity distribution from the prescribed pixel
                displayvtvf(vt,optr,maxveloc,dbrange,X,Y,slice,displayvfspace)
            end
        end
    end
end