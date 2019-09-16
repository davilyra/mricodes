% Reconstruction of the time-velocity distribuition from multiple coils
% and defined points in a voxel in a slice, using sos nuFFT and SPIRiT
% algorithm. The points correspond to: Left Carotid Bifurcation, Left
% Jugular Vein, Right External Carotid Artery, Right Internal Carotid
% Artery and Right Jugular Vein.
%
% By Davi Marco Lyra-Leite <davi@ieee.org>
% November 5th, 2011
%
% Edited by Joao Carvalho on Nov 6th, 2011.
% Edited by Davi Lyra-Leite on March 11th, 2012, and October 26th, 2012.


% cleaning previous data
clear,clc,close all

%% loading data
%slice = input('Com qual slice você deseja trabalhar (1 a 5)? \n R: ','s');
slice = '4';

% selecting the viewordering method
viewordering = '0';

%% analysis points (for slice 4)
axesX = round([74.8667 77.9333 35 33.4667 28.8667 41.6444 65.6667]);
axesY = round([71.0333 66.9444 75.1222 70.5222 65.9222 58.2556 58.2556]);

%% undersampling factors
usfactors = [1 2 4];
% usfactors = [1 4];


%% initializing data matrices

% SPIRiT
xyspi = zeros(115,115,3);
vtspi = zeros(32,43,7,3);
vfspi = zeros(32,43,7,3);
xyerrorspi = zeros(115,115,3);
vterrorspi = zeros(32,43,7,3);
vferrorspi = zeros(32,43,7,3);

% SOS
xysos = zeros(115,115,3);
vtsos = zeros(32,43,7,3);
vfsos = zeros(32,43,7,3);
xyerrorsos = zeros(115,115,3);
vterrorsos = zeros(32,43,7,3);
vferrorsos = zeros(32,43,7,3);


%% loading data
disp('loading data and calculating errors...')
circmask = mask(115); 
% circmask  = ones(115,115);
cd datapaper/
% for usf = 1:2,%
for usf = 1:3
    usfactorstr = num2str(usfactors(usf));
    
    disp(sprintf('under sampling factor: %d',usfactors(usf)));pause(.1)

    % loading SOS data
    if usf == 1
        cd sosnufft_3 % loading the original fully-sampled reference
        load(['slice_',slice,'_sosnufft_',usfactorstr,'.mat'])
%         xy = abs(xy); xyvt = abs(xyvt);
%         xy = (xy)/max(xy(:)); xyvt = (xyvt)/max(xyvt(:));
        xysos(:,:,usf) = xy.*circmask;
        for voxel=1:size(axesX,2),
            x = axesX(voxel); y = axesY(voxel);
            vtsos(:,:,voxel,usf) = flipud(permute(abs(xyvt(x,y,:,:)),[3 4 1 2]));
            vfsos(:,:,voxel,usf) = abs(fftshift(abs(fft(flipud(permute(xyvt(x,y,:,:),[3 4 1 2])),[],2)),2));
        end;
    else
        cd spirit_JL1
        load(['recon_slice_sos_slice_',slice,'_us_',usfactorstr,'_method_',viewordering,'.mat'])
%         xy = abs(xy); xyvt = abs(xyvt);
%         xy = (xy)/max(xy(:)); xyvt = (xyvt)/max(xyvt(:));
        xysos(:,:,usf) = xy.*circmask;
        for voxel=1:size(axesX,2),
            x = axesX(voxel); y = axesY(voxel);
            vtsos(:,:,voxel,usf) = flipud(permute(abs(xyvt(x,y,:,:)),[3 4 1 2]));
            vfsos(:,:,voxel,usf) = abs(fftshift(abs(fft(flipud(permute(xyvt(x,y,:,:),[3 4 1 2])),[],2)),2));
        end;
    end
    
    
    if usf ~=1
        % loading SPIRiT data
        load(['recon_spirit_slice_',slice,'_us_',usfactorstr,'_method_',viewordering,'.mat'])
%         xy = abs(xy); xyvt = abs(xyvt);
%         xy = (xy)/max(xy(:)); xyvt = (xyvt)/max(xyvt(:));
        xyspi(:,:,usf) = xy.*circmask;
        for voxel=1:size(axesX,2),
            x = axesX(voxel); y = axesY(voxel);
            vtspi(:,:,voxel,usf) = flipud(permute(abs(xyvt(x,y,:,:)),[3 4 1 2]));
            vfspi(:,:,voxel,usf) = abs(fftshift(abs(fft(flipud(permute(xyvt(x,y,:,:),[3 4 1 2])),[],2)),2));
        end;
    end
    cd ..
    vfsos(:,:,:,usf) = vfsos(:,:,:,usf)/max(max(max(vfsos(:,:,:,usf))));
    vfspi(:,:,:,usf) = vfspi(:,:,:,usf)/max(max(max(vfspi(:,:,:,usf))));
    
    % reference images (fully-sampled SOS)
    if usf == 1,        
        xyref = xysos(:,:,1);
        vtref = vtsos(:,:,:,1);
        vfref = vfsos(:,:,:,1);
    end;
      
    % calculating errors
    xyerrorsos(:,:,usf) = xyref - xysos(:,:,usf);
    xyerrorspi(:,:,usf) = xyref - xyspi(:,:,usf);
    for voxel = 1:size(axesX,2),
        % vt error
        vterrorsos(:,:,voxel,usf) = vtref(:,:,voxel) - vtsos(:,:,voxel,usf);
        vterrorspi(:,:,voxel,usf) = vtref(:,:,voxel) - vtspi(:,:,voxel,usf);
        
        % vferror
        vferrorsos(:,:,voxel,usf) = vfref(:,:,voxel) - vfsos(:,:,voxel,usf);
        vferrorspi(:,:,voxel,usf) = vfref(:,:,voxel) - vfspi(:,:,voxel,usf);
    end;
end;
clear xyvt xy circmask
cd ..


%% calculating the SER
disp('calculating the signal-to-error ratios...');pause(.1)
serdb1 = zeros(16,3);
% for usf = 1:2,%
for usf = 1:3
    serdb1(1,usf) = serdb(xyref,xysos(:,:,usf));
    serdb1(2,usf) = serdb(xyref,xyspi(:,:,usf));
    for voxel = 1:size(axesX,2),
        serdb1(voxel+2,usf) = serdb(vtref(:,:,voxel),vtsos(:,:,voxel,usf));
        serdb1(voxel+9,usf) = serdb(vtref(:,:,voxel),vtspi(:,:,voxel,usf));
    end;
end;
cd ser
save(['ser_spirit3D_slice_',slice,'_method_',viewordering,'.mat'], 'serdb1');
cd ..
serdb1,


%% normalizing images before writing to file
disp('preparing images for writing...')
maxxy = max(abs(xyref(:)));
xysos = abs(xysos)/maxxy;
xyspi = abs(xyspi)/maxxy;
xyerrorsos = abs(xyerrorsos)/maxxy;
xyerrorspi = abs(xyerrorspi)/maxxy;
for voxel = 1 : size(axesX,2),
    maxvt = max(max(abs(vtref(:,:,voxel))));
    vtsos(:,:,voxel,:) = abs(vtsos(:,:,voxel,:))/maxvt;
    vtspi(:,:,voxel,:) = abs(vtspi(:,:,voxel,:))/maxvt;
    vterrorsos(:,:,voxel,:) = abs(vterrorsos(:,:,voxel,:))/maxvt;
    vterrorspi(:,:,voxel,:) = abs(vterrorspi(:,:,voxel,:))/maxvt;
    
    maxvf = max(max(abs(vfref(:,:,voxel))));
    vfsos(:,:,voxel,:) = abs(vfsos(:,:,voxel,:))/maxvf;
    vfspi(:,:,voxel,:) = abs(vfspi(:,:,voxel,:))/maxvf;
    vferrorsos(:,:,voxel,:) = abs(vferrorsos(:,:,voxel,:))/maxvf;
    vferrorspi(:,:,voxel,:) = abs(vferrorspi(:,:,voxel,:))/maxvf;
end;


%% saving images
disp('saving images...')
cd images/spirit3D_vf % inversion data using spirit and nufft
% for usf = 1:2,%
for usf = 1:3
    usfactorstr = num2str(usfactors(usf));
    
%     % saving xy images
    imwrite(xysos(:,:,usf),['xysos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'.tiff'],'TIFF');
%     
%     % saving xy error images
    imwrite(xyerrorsos(:,:,usf),['xyerrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'.tiff'],'TIFF');
%     
%     % saving vt images (SOS)
    imwrite(vtsos(:,:,1,usf),['vtsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LCB.tiff'],'TIFF');
    imwrite(vtsos(:,:,2,usf),['vtsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LJV.tiff'],'TIFF');
    imwrite(vtsos(:,:,3,usf),['vtsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RECA.tiff'],'TIFF');
    imwrite(vtsos(:,:,4,usf),['vtsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RICA.tiff'],'TIFF');
    imwrite(vtsos(:,:,5,usf),['vtsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RJV.tiff'],'TIFF');
    imwrite(vtsos(:,:,6,usf),['vtsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RVA.tiff'],'TIFF');
    imwrite(vtsos(:,:,7,usf),['vtsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LVA.tiff'],'TIFF');
%     
%     % saving vt error images (SOS)
    imwrite(vterrorsos(:,:,1,usf),['vterrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LCB.tiff'],'TIFF');
    imwrite(vterrorsos(:,:,2,usf),['vterrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LJV.tiff'],'TIFF');
    imwrite(vterrorsos(:,:,3,usf),['vterrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RECA.tiff'],'TIFF');
    imwrite(vterrorsos(:,:,4,usf),['vterrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RICA.tiff'],'TIFF');
    imwrite(vterrorsos(:,:,5,usf),['vterrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RJV.tiff'],'TIFF');
    imwrite(vterrorsos(:,:,6,usf),['vterrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RVA.tiff'],'TIFF');
    imwrite(vterrorsos(:,:,7,usf),['vterrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LVA.tiff'],'TIFF');
    

%     % saving vf images (SOS)
%     imwrite(vfsos(:,:,1,usf),['vfsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LCB.tiff'],'TIFF');
%     imwrite(vfsos(:,:,2,usf),['vfsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LJV.tiff'],'TIFF');
%     imwrite(vfsos(:,:,3,usf),['vfsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RECA.tiff'],'TIFF');
%     imwrite(vfsos(:,:,4,usf),['vfsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RICA.tiff'],'TIFF');
%     imwrite(vfsos(:,:,5,usf),['vfsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RJV.tiff'],'TIFF');
%     imwrite(vfsos(:,:,6,usf),['vfsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RVA.tiff'],'TIFF');
%     imwrite(vfsos(:,:,7,usf),['vfsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LVA.tiff'],'TIFF');
%     
%         
% 	if usf ~= 1
%         % saving vf error images (SOS)
%         imwrite(vferrorsos(:,:,1,usf),['vferrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LCB.tiff'],'TIFF');
%         imwrite(vferrorsos(:,:,2,usf),['vferrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LJV.tiff'],'TIFF');
%         imwrite(vferrorsos(:,:,3,usf),['vferrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RECA.tiff'],'TIFF');
%         imwrite(vferrorsos(:,:,4,usf),['vferrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RICA.tiff'],'TIFF');
%         imwrite(vferrorsos(:,:,5,usf),['vferrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RJV.tiff'],'TIFF');
%         imwrite(vferrorsos(:,:,6,usf),['vferrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RVA.tiff'],'TIFF');
%         imwrite(vferrorsos(:,:,7,usf),['vferrorsos_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LVA.tiff'],'TIFF');
% 	end

    if usf ~= 1
%         % saving xy images
        imwrite(xyspi(:,:,usf),['xyspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'.tiff'],'TIFF');
%         
%         % saving xy error images
        imwrite(xyerrorspi(:,:,usf),['xyerrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'.tiff'],'TIFF');
%         
%         % saving vt images (SPIRiT)
        imwrite(vtspi(:,:,1,usf),['vtspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LCB.tiff'],'TIFF');
        imwrite(vtspi(:,:,2,usf),['vtspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LJV.tiff'],'TIFF');
        imwrite(vtspi(:,:,3,usf),['vtspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RECA.tiff'],'TIFF');
        imwrite(vtspi(:,:,4,usf),['vtspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RICA.tiff'],'TIFF');
        imwrite(vtspi(:,:,5,usf),['vtspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RJV.tiff'],'TIFF');
        imwrite(vtspi(:,:,6,usf),['vtspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RVA.tiff'],'TIFF');
        imwrite(vtspi(:,:,7,usf),['vtspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LVA.tiff'],'TIFF');
%         
%         % saving vt error images (SPIRiT)
        imwrite(vterrorspi(:,:,1,usf),['vterrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LCB.tiff'],'TIFF');
        imwrite(vterrorspi(:,:,2,usf),['vterrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LJV.tiff'],'TIFF');
        imwrite(vterrorspi(:,:,3,usf),['vterrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RECA.tiff'],'TIFF');
        imwrite(vterrorspi(:,:,4,usf),['vterrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RICA.tiff'],'TIFF');
        imwrite(vterrorspi(:,:,5,usf),['vterrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RJV.tiff'],'TIFF');
        imwrite(vterrorspi(:,:,6,usf),['vterrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RVA.tiff'],'TIFF');
        imwrite(vterrorspi(:,:,7,usf),['vterrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LVA.tiff'],'TIFF');
        
%         % saving vf images (SPIRiT)
%         imwrite(vfspi(:,:,1,usf),['vfspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LCB.tiff'],'TIFF');
%         imwrite(vfspi(:,:,2,usf),['vfspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LJV.tiff'],'TIFF');
%         imwrite(vfspi(:,:,3,usf),['vfspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RECA.tiff'],'TIFF');
%         imwrite(vfspi(:,:,4,usf),['vfspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RICA.tiff'],'TIFF');
%         imwrite(vfspi(:,:,5,usf),['vfspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RJV.tiff'],'TIFF');
%         imwrite(vfspi(:,:,6,usf),['vfspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RVA.tiff'],'TIFF');
%         imwrite(vfspi(:,:,7,usf),['vfspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LVA.tiff'],'TIFF');
%         
%         % saving vf error images (SPIRiT)
%         imwrite(vferrorspi(:,:,1,usf),['vferrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LCB.tiff'],'TIFF');
%         imwrite(vferrorspi(:,:,2,usf),['vferrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LJV.tiff'],'TIFF');
%         imwrite(vferrorspi(:,:,3,usf),['vferrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RECA.tiff'],'TIFF');
%         imwrite(vferrorspi(:,:,4,usf),['vferrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RICA.tiff'],'TIFF');
%         imwrite(vferrorspi(:,:,5,usf),['vferrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RJV.tiff'],'TIFF');
%         imwrite(vferrorspi(:,:,6,usf),['vferrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_RVA.tiff'],'TIFF');
%         imwrite(vferrorspi(:,:,7,usf),['vferrorspi_sl_',slice,'_us_',usfactorstr,'_method_',viewordering,'_LVA.tiff'],'TIFF');
    end
end;
cd ..
cd ..

disp('done!')