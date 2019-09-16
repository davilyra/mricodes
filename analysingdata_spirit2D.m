% Reconstruction of the time-velocity distribuition from multiple coils
% and defined points in a voxel in a slice, using sos nuFFT and SPIRiT
% algorithm. The points correspond to: Left Carotid Bifurcation, Left
% Jugular Vein, Right External Carotid Artery, Right Internal Carotid
% Artery and Right Jugular Vein.
%
% By Davi Marco Lyra-Leite <davi@ieee.org>
% November 5th, 2011 and February 14th
%
% Edited by Joao Carvalho on Nov 6th, 2011.
% Edited by Davi Lyra-Leite Carvalho on March 11th, 2012.
%
% In order to perform studies evaluating data without spiral FVE steps
% inversion use data from: spirit_2 and sosnufft_3


% cleaning previous data
clear,clc,close all

%% loading data
%slice = input('Com qual slice você deseja trabalhar (1 a 5)? \n R: ','s');
slice = '4';

%% analysis points
axesX = round([74.8667 77.9333 35 33.4667 28.8667 41.6444 65.6667]);
axesY = round([71.0333 66.9444 75.1222 70.5222 65.9222 58.2556 58.2556]);

%% undersampling factors
usfactors = [1 2 4];

%% initializing data matrices

% SPIRiT
xyspi = zeros(115,115,3);
vtspi = zeros(32,43,7,3);
xyerrorspi = zeros(115,115,3);
vterrorspi = zeros(32,43,7,3);

% SOS
xysos = zeros(115,115,3);
vtsos = zeros(32,43,7,3);
xyerrorsos = zeros(115,115,3);
vterrorsos = zeros(32,43,7,3);

%% loading data
disp('loading data and calculating errors...')
circmask = mask(115); 
% circmask  = ones(115,115);
cd datapaper\
for usf = 1:3,
    usfactorstr = num2str(usfactors(usf));
    
    disp(sprintf('under sampling factor: %d',usfactors(usf)));pause(.1)

    accelfactor = double(usfactors(usf));
    scalefactor= sqrt(accelfactor);
    
    % loading SOS data
%     cd sosnufft_2 % sos without inversion
    cd sosnufft_3 % sos with inversion
%     cd sosnufft_coil_selection

    load(['slice_',slice,'_sosnufft_',usfactorstr,'.mat'])
%     load(['slice_',slice,'_sosnufft_coils_23_',usfactorstr,'.mat'])
%     load(['slice_',slice,'_sosnufft_coils_12_',usfactorstr,'.mat'])
%     xy = accelfactor*abs(xy); xyvt = accelfactor*abs(xyvt);
    xy = abs(xy); xyvt = abs(xyvt);
    xysos(:,:,usf) = xy.*circmask;
    for voxel=1:size(axesX,2),
        x = axesX(voxel); y = axesY(voxel);
        vtsos(:,:,voxel,usf) = flipud(permute(xyvt(x,y,:,:),[3 4 1 2]));
    end;
    cd ..
    
    
    % loading SPIRiT data
%     cd spirit % spirit without inversion
%     cd spirit_2 % spirit with inversion
%     cd spirit_coil_selection
    cd spirit_JL1

%     load(['slice_',slice,'_',usfactorstr,'.mat']) % for spirit data with/without inversion
%     load(['slice_',slice,'_',usfactorstr,'_xy.mat']) % for spirita data without inversion
    load(['recon_slice_',slice,'_',usfactorstr,'.mat']) % for spirita data without inversion
%     load(['slice_',slice,'_spirit_coils_23_',usfactorstr,'.mat'])
%     load(['slice_',slice,'_spirit_coils_12_',usfactorstr,'.mat'])
    xy = abs(xy); xyvt = abs(xyvt);
    xyspi(:,:,usf) = xy.*circmask;
    for voxel=1:size(axesX,2),
        x = axesX(voxel); y = axesY(voxel);
        vtspi(:,:,voxel,usf) = flipud(permute(xyvt(x,y,:,:),[3 4 1 2]));
    end;
    cd ..
    
    % reference images (fully-sampled SOS)
    if usf==1,        
        xyref = xysos(:,:,1);
        vtref = vtsos(:,:,:,1);
    end;
      
    %calculating errors
    xyerrorsos(:,:,usf) = xyref - xysos(:,:,usf);
    xyerrorspi(:,:,usf) = xyref - xyspi(:,:,usf);
    for voxel=1:size(axesX,2),
        vterrorsos(:,:,voxel,usf) = vtref(:,:,voxel) - vtsos(:,:,voxel,usf);
        vterrorspi(:,:,voxel,usf) = vtref(:,:,voxel) - vtspi(:,:,voxel,usf);
    end;
end;
clear xyvt xy circmask
cd ..


%% calculating the SER
disp('calculating the signal-to-error ratios...');pause(.1)
ser = zeros(16,3);
for usf = 1:3,
    ser(1,usf) = serdb(xyref,xysos(:,:,usf));
    ser(2,usf) = serdb(xyref,xyspi(:,:,usf));
    for voxel = 1:size(axesX,2),
        ser(voxel+2,usf) = serdb(vtref(:,:,voxel),vtsos(:,:,voxel,usf));
        ser(voxel+9,usf) = serdb(vtref(:,:,voxel),vtspi(:,:,voxel,usf));
    end;
end;
cd ser
% save(['ser_spirit_slice_',slice,'.mat'], 'ser');
% save(['ser_spirit_2_slice_',slice,'.mat'], 'ser');
% save(['ser_spirit_coils23_slice_',slice,'.mat'], 'ser');
% save(['ser_spirit_coils12_slice_',slice,'.mat'], 'ser');
save(['ser_spirit_JL1_slice_',slice,'.mat'], 'ser');
cd ..
ser,


% %% plotting couple of results
% figure
% subplot(236), imagesc(abs(xyerrorspi(:,:,3))), colormap(gray), colorbar, axis off
% title('error spirit 4-fold');
% subplot(235), imagesc(abs(xyerrorspi(:,:,2))), colormap(gray), colorbar, axis off
% title('error spirit 2-fold');
% subplot(234), imagesc(abs(xyerrorspi(:,:,1))), colormap(gray), colorbar, axis off
% title('error spirit fully-sampled');
% subplot(233), imagesc(abs(xyerrorsos(:,:,3))), colormap(gray), colorbar, axis off
% title('error sos 4-fold');
% subplot(232), imagesc(abs(xyerrorsos(:,:,2))), colormap(gray), colorbar, axis off
% title('error sos 2-fold');
% subplot(231), imagesc(abs(xysos(:,:,1))), colormap(gray), colorbar, axis off
% title('sos fully sampled');


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
end;


%% saving images
disp('saving images...')
% cd images\spirit % data using spirit and nufft
% cd images\spirit_inversion % inversion data using spirit and nufft
% cd images\spirit_coil_selection % inversion data using only coils 2 and 3 for spirit and nufft
% cd images\spirit_coil_selection_2 % inversion data using only coils 2 and 3 for spirit and nufft
cd images\spirit_JL1 % inversion data using spirit and nufft
for usf = 1:3,
    usfactorstr = num2str(usfactors(usf));
    
    % saving xy images
    imwrite(xysos(:,:,usf),['xysos_sl',slice,'_us',usfactorstr,'.png'],'PNG');
    imwrite(xyspi(:,:,usf),['xyspi_sl',slice,'_us',usfactorstr,'.png'],'PNG');

    % saving vt images (SOS)
    imwrite(vtsos(:,:,1,usf),['vtsos_sl',slice,'_us',usfactorstr,'_LCB.png'],'PNG');
    imwrite(vtsos(:,:,2,usf),['vtsos_sl',slice,'_us',usfactorstr,'_LJV.png'],'PNG');
    imwrite(vtsos(:,:,3,usf),['vtsos_sl',slice,'_us',usfactorstr,'_RECA.png'],'PNG');
    imwrite(vtsos(:,:,4,usf),['vtsos_sl',slice,'_us',usfactorstr,'_RICA.png'],'PNG');
    imwrite(vtsos(:,:,5,usf),['vtsos_sl',slice,'_us',usfactorstr,'_RJV.png'],'PNG');
    imwrite(vtsos(:,:,6,usf),['vtsos_sl',slice,'_us',usfactorstr,'_RVA.png'],'PNG');
    imwrite(vtsos(:,:,7,usf),['vtsos_sl',slice,'_us',usfactorstr,'_LVA.png'],'PNG');
    
    % saving vt images (SPIRiT)
    imwrite(vtspi(:,:,1,usf),['vtspi_sl',slice,'_us',usfactorstr,'_LCB.png'],'PNG');
    imwrite(vtspi(:,:,2,usf),['vtspi_sl',slice,'_us',usfactorstr,'_LJV.png'],'PNG');
    imwrite(vtspi(:,:,3,usf),['vtspi_sl',slice,'_us',usfactorstr,'_RECA.png'],'PNG');
    imwrite(vtspi(:,:,4,usf),['vtspi_sl',slice,'_us',usfactorstr,'_RICA.png'],'PNG');
    imwrite(vtspi(:,:,5,usf),['vtspi_sl',slice,'_us',usfactorstr,'_RJV.png'],'PNG');
    imwrite(vtspi(:,:,6,usf),['vtspi_sl',slice,'_us',usfactorstr,'_RVA.png'],'PNG');
    imwrite(vtspi(:,:,7,usf),['vtspi_sl',slice,'_us',usfactorstr,'_LVA.png'],'PNG');
    
    % saving xy error images
    imwrite(xyerrorsos(:,:,usf),['xyerrorsos_sl',slice,'_us',usfactorstr,'.png'],'PNG');
    imwrite(xyerrorspi(:,:,usf),['xyerrorspi_sl',slice,'_us',usfactorstr,'.png'],'PNG');

    % saving vt error images (SOS)
    imwrite(vterrorsos(:,:,1,usf),['vterrorsos_sl',slice,'_us',usfactorstr,'_LCB.png'],'PNG');
    imwrite(vterrorsos(:,:,2,usf),['vterrorsos_sl',slice,'_us',usfactorstr,'_LJV.png'],'PNG');
    imwrite(vterrorsos(:,:,3,usf),['vterrorsos_sl',slice,'_us',usfactorstr,'_RECA.png'],'PNG');
    imwrite(vterrorsos(:,:,4,usf),['vterrorsos_sl',slice,'_us',usfactorstr,'_RICA.png'],'PNG');
    imwrite(vterrorsos(:,:,5,usf),['vterrorsos_sl',slice,'_us',usfactorstr,'_RJV.png'],'PNG');
    imwrite(vterrorsos(:,:,6,usf),['vterrorsos_sl',slice,'_us',usfactorstr,'_RVA.png'],'PNG');
    imwrite(vterrorsos(:,:,7,usf),['vterrorsos_sl',slice,'_us',usfactorstr,'_LVA.png'],'PNG');
    
    % saving vt error images (SPIRiT)
    imwrite(vterrorspi(:,:,1,usf),['vterrorspi_sl',slice,'_us',usfactorstr,'_LCB.png'],'PNG');
    imwrite(vterrorspi(:,:,2,usf),['vterrorspi_sl',slice,'_us',usfactorstr,'_LJV.png'],'PNG');
    imwrite(vterrorspi(:,:,3,usf),['vterrorspi_sl',slice,'_us',usfactorstr,'_RECA.png'],'PNG');
    imwrite(vterrorspi(:,:,4,usf),['vterrorspi_sl',slice,'_us',usfactorstr,'_RICA.png'],'PNG');
    imwrite(vterrorspi(:,:,5,usf),['vterrorspi_sl',slice,'_us',usfactorstr,'_RJV.png'],'PNG');
    imwrite(vterrorspi(:,:,6,usf),['vterrorspi_sl',slice,'_us',usfactorstr,'_RVA.png'],'PNG');
    imwrite(vterrorspi(:,:,7,usf),['vterrorspi_sl',slice,'_us',usfactorstr,'_LVA.png'],'PNG');
end;
cd ..
cd ..

disp('done!')