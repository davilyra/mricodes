% Reconstruction of the time-velocity distribuition from multiple coils
% and defined points in a voxel in a slice, using sos nuFFT and SPIRiT
% algorithm. The points correspond to: Left Carotid Bifurcation, Left
% Jugular Vein, Right External Carotid Artery, Right Internal Carotid
% Artery and Right Jugular Vein.
%
% By Davi Marco Lyra Leite <davi@ieee.org>
% November 5, 2011

% cleaning previous data
clear,clc,close all


%% analysis points
axesX = [74.8667 77.9333 35 33.4667 28.8667];
axesY = [71.0333 66.9444 75.1222 70.5222 65.9222];


%% loading data
t = input('Com qual slice você deseja trabalhar (1 a 5)? \n R: ','s');

%% initializing data matrices
% SPIRiT
xyvtspi = zeros(115,115,32,43,3);
xyspi = zeros(115,115,3);
xyspisave = zeros(115,115,3);
xyspien = zeros(115,115,3);
xyspienc = zeros(115,115,3);
vtspi = zeros(32,43,5,3);

% nuFFT SOS
xyvt = zeros(115,115,32,43,3);
xy = zeros(115,115,3);
xysave = zeros(115,115,3);
xyen = zeros(115,115,3);
xyenc = zeros(115,115,3);
vt = zeros(32,43,5,3);

% error
error0 = zeros(115,115,3);
error1 = zeros(115,115,3);

% flow profiles
vt_en = zeros(32,43,5,3);
vtspi_en = zeros(32,43,5,3);
error_vt = zeros(32,43,5,3);
error_vtsos = zeros(32,43,5,3);

% circular mask
circmask = mask(115);

% undersampling factors
t1 = [1 2 4];

for i = 1:3
    ta1 = num2str(t1(i));

    disp('loading data...')
    cd datapaper\
    % for SPIRiT data
    load(['slice_',t,'_',ta1,'.mat'])
    load(['slice_',t,'_',ta1,'_xy.mat'])

    xyvtspi(:,:,:,:,i) = xyvt;
    xyspi(:,:,i) = xy;

    % for nuFFT SOS data
    load(['slice_',t,'_sosnufft_',ta1,'.mat'])
    xyvt(:,:,:,:,i) = xyvt;
    xy(:,:,i) = xy;
    
    cd ..
    
    %% takes the velocity distribution from the prescribed pixel
    % SPIRiT
    vtspi_data = xyvtspi(round(axesX(1)),round(axesY(1)),:,:,i); % Left Carotid Bifurcation
    vtspi(:,:,1,i) = permute(vtspi_data,[3 4 1 2]);

    vtspi_data = xyvtspi(round(axesX(2)),round(axesY(2)),:,:,i); % Left Jugular Vein
    vtspi(:,:,2,i) = permute(vtspi_data,[3 4 1 2]);

    vtspi_data = xyvtspi(round(axesX(3)),round(axesY(3)),:,:,i); % Right External Carotid Artery
    vtspi(:,:,3,i) = permute(vtspi_data,[3 4 1 2]);

    vtspi_data = xyvtspi(round(axesX(4)),round(axesY(4)),:,:,i); % Right Internal Carotid Artery
    vtspi(:,:,4,i) = permute(vtspi_data,[3 4 1 2]);

    vtspi_data = xyvtspi(round(axesX(5)),round(axesY(5)),:,:,i); % Right Jugular Vein
    vtspi(:,:,5,i) = permute(vtspi_data,[3 4 1 2]);

    clear vtspi_data

    
    % nuFFT SOS
    vt_data = xyvt(round(axesX(1)),round(axesY(1)),:,:,i); % Left Carotid Bifurcation
    vt(:,:,1,i) = permute(vt_data,[3 4 1 2]);

    vt_data = xyvt(round(axesX(2)),round(axesY(2)),:,:,i); % Left Jugular Vein
    vt(:,:,2,i) = permute(vt_data,[3 4 1 2]);

    vt_data = xyvt(round(axesX(3)),round(axesY(3)),:,:,i); % Right External Carotid Artery
    vt(:,:,3,i) = permute(vt_data,[3 4 1 2]);

    vt_data = xyvt(round(axesX(4)),round(axesY(4)),:,:,i); % Right Internal Carotid Artery
    vt(:,:,4,i) = permute(vt_data,[3 4 1 2]);

    vt_data = xyvt(round(axesX(5)),round(axesY(5)),:,:,i); % Right Jugular Vein
    vt(:,:,5,i) = permute(vt_data,[3 4 1 2]);

    clear vt_data


    %% errors
    disp('calculating erros...')
        
    % SPIRiT and nuFFT
    xyspien(:,:,i) = (xyspi(:,:,i))/energyd(abs(xyspi(:,:,i)));
    xyen(:,:,i) = (xy(:,:,i))/energyd(abs(xy(:,:,i)));
    
    xyenc(:,:,i) = xyen(:,:,i).*circmask;
    xyspienc(:,:,i) = xyen(:,:,i).*circmask;
    
    error0(:,:,i) = abs(xyenc(:,:,1) - xyspienc(:,:,i));
    error1(:,:,i) = abs(xyenc(:,:,1) - xyenc(:,:,i));

    for j = 1 : 5
        vt_en(:,:,j,i) = (vt(:,:,j,i))/energyd(abs(vt(:,:,j,i)));
        vtspi_en(:,:,j,i) = (vtspi(:,:,j,i))/energyd(abs(vtspi(:,:,j,i)));
        error_vt(:,:,j,i) = abs(vt_en(:,:,j,1) - vtspi_en(:,:,j,i));
        error_vtsos(:,:,j,i) = abs(vt_en(:,:,j,1) - vt_en(:,:,j,i));
    end
    
    
    %%  saving into images
    cd images\
    disp('saving images...')
    disp(' ')
    
    % saving data (SPIRiT)
    xyspisave(:,:,i) = abs(xyspi(:,:,i))/(max(max(abs(xyspi(:,:,i)))));
    imwrite(abs(xyspisave(:,:,i)),['xyspi_',t,'_accel_',ta1,'.png'],'PNG');

    % saving data (nuFFT SOS)
    xysave(:,:,i) = abs(xy(:,:,i))/max(max(abs(xy(:,:,i))));
    imwrite(abs(xysave(:,:,i)),['xynufft_',t,'_accel_',ta1,'.png'],'PNG');

    % SPIRiT
    imwrite((flipud(abs(vtspi(:,:,1,i))/max(max(abs(vtspi(:,:,1,i)))))),['vtspi_',t,'_',ta1,'_LCB.png'],'PNG');
    imwrite((flipud(abs(vtspi(:,:,2,i))/max(max(abs(vtspi(:,:,2,i)))))),['vtspi_',t,'_',ta1,'_LJV.png'],'PNG');
    imwrite((flipud(abs(vtspi(:,:,3,i))/max(max(abs(vtspi(:,:,3,i)))))),['vtspi_',t,'_',ta1,'_RECA.png'],'PNG');
    imwrite((flipud(abs(vtspi(:,:,4,i))/max(max(abs(vtspi(:,:,4,i)))))),['vtspi_',t,'_',ta1,'_RICA.png'],'PNG');
    imwrite((flipud(abs(vtspi(:,:,5,i))/max(max(abs(vtspi(:,:,5,i)))))),['vtspi_',t,'_',ta1,'_RJV.png'],'PNG');

    % nuFFT SOS
    imwrite((flipud(abs(vt(:,:,1,i))/max(max(abs(vt(:,:,1,i)))))),['vtnufft_',t,'_',ta1,'_LCB.png'],'PNG');
    imwrite((flipud(abs(vt(:,:,2,i))/max(max(abs(vt(:,:,2,i)))))),['vtnufft_',t,'_',ta1,'_LJV.png'],'PNG');
    imwrite((flipud(abs(vt(:,:,3,i))/max(max(abs(vt(:,:,3,i)))))),['vtnufft_',t,'_',ta1,'_RECA.png'],'PNG');
    imwrite((flipud(abs(vt(:,:,4,i))/max(max(abs(vt(:,:,4,i)))))),['vtnufft_',t,'_',ta1,'_RICA.png'],'PNG');
    imwrite((flipud(abs(vt(:,:,5,i))/max(max(abs(vt(:,:,5,i)))))),['vtnufft_',t,'_',ta1,'_RJV.png'],'PNG');

    % errors
    imwrite((flipud(abs(error0(:,:,i))/max(abs(error0(:))))),['error_nufftfs_spirit_',t,'_us_',ta1,'.png'],'PNG');
    imwrite((flipud(abs(error1(:,:,i))/max(abs(error1(:))))),['error_nufftfs_nufft_',t,'_us_',ta1,'.png'],'PNG');
    imwrite((flipud(abs(error_vt(:,:,1,i))/max(max(abs(error_vt(:,:,1,i)))))),['error_spirit_nufft_',t,'_us_',ta1,'_LCB.png'],'PNG');
    imwrite((flipud(abs(error_vt(:,:,2,i))/max(max(abs(error_vt(:,:,2,i)))))),['error_spirit_nufft_',t,'_us_',ta1,'_LJV.png'],'PNG');
    imwrite((flipud(abs(error_vt(:,:,3,i))/max(max(abs(error_vt(:,:,3,i)))))),['error_spirit_nufft_',t,'_us_',ta1,'_RECA.png'],'PNG');
    imwrite((flipud(abs(error_vt(:,:,4,i))/max(max(abs(error_vt(:,:,4,i)))))),['error_spirit_nufft_',t,'_us_',ta1,'_RICA.png'],'PNG');
    imwrite((flipud(abs(error_vt(:,:,5,i))/max(max(abs(error_vt(:,:,5,i)))))),['error_spirit_nufft_',t,'_us_',ta1,'_RJV.png'],'PNG');
    
    imwrite((flipud(abs(error_vtsos(:,:,1,i))/max(max(abs(error_vtsos(:,:,1,i)))))),['error_nufft_nufft_',t,'_us_',ta1,'_LCB.png'],'PNG');
    imwrite((flipud(abs(error_vtsos(:,:,2,i))/max(max(abs(error_vtsos(:,:,2,i)))))),['error_nufft_nufft_',t,'_us_',ta1,'_LJV.png'],'PNG');
    imwrite((flipud(abs(error_vtsos(:,:,3,i))/max(max(abs(error_vtsos(:,:,3,i)))))),['error_nufft_nufft_',t,'_us_',ta1,'_RECA.png'],'PNG');
    imwrite((flipud(abs(error_vtsos(:,:,4,i))/max(max(abs(error_vtsos(:,:,4,i)))))),['error_nufft_nufft_',t,'_us_',ta1,'_RICA.png'],'PNG');
    imwrite((flipud(abs(error_vtsos(:,:,5,i))/max(max(abs(error_vtsos(:,:,5,i)))))),['error_nufft_nufft_',t,'_us_',ta1,'_RJV.png'],'PNG');
    
    cd ..

end

%% calculating the SER
disp('calculating SERs...')
ser = zeros(12,3);
ser(1,1) = snr(xyenc(:,:,1),xyenc(:,:,1));
ser(1,2) = snr(xyenc(:,:,1),xyenc(:,:,2));
ser(1,3) = snr(xyenc(:,:,1),xyenc(:,:,3));
ser(2,1) = snr(xyenc(:,:,1),xyspienc(:,:,1));
ser(2,2) = snr(xyenc(:,:,1),xyspienc(:,:,2));
ser(2,3) = snr(xyenc(:,:,1),xyspienc(:,:,3));

for j = 3:7
    ser(j,1) = snr(vt_en(:,:,(j - 2),1),vtspi_en(:,:,(j - 2),1));
    ser(j,2) = snr(vt_en(:,:,(j - 2),1),vtspi_en(:,:,(j - 2),2));
    ser(j,3) = snr(vt_en(:,:,(j - 2),1),vtspi_en(:,:,(j - 2),3));
end

for j = 8:12
    ser(j,1) = snr(vt_en(:,:,(j - 7),1),vt_en(:,:,(j - 7),1));
    ser(j,2) = snr(vt_en(:,:,(j - 7),1),vt_en(:,:,(j - 7),2));
    ser(j,3) = snr(vt_en(:,:,(j - 7),1),vt_en(:,:,(j - 7),3));
end

save(['ser_slice_',t,'.mat'], 'ser');
