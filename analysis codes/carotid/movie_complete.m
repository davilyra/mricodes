%% Clearing the Workspace
clear all; close all; clc

accelfactor = 8; % CS reduction factor
frames = [1:9 9];
nframes = length(frames); % number of BBTI's

folder2 = [400 500 600 700 800 900 300 1000 1100]; % BBTI values from the performed study
folder1 = [10 11 12 13 14 15 16 17 18]; % Datasets id
folder = [folder1' folder2'];
folder = sortrows(folder,2);
clear folder1 folder2

load('/home/mri-vh/Davi/carotid/params.mat')

movie1 = moviein(nframes);

for ii = 1:nframes;
	load(['accelfactor_',num2str(accelfactor),'/carotid_fullysampled.mat']);
    mip = max(squeeze(rd_full(:,:,:,frames(ii))), [], 3);
    subplot(231)
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(frames(ii),2)),' ms'])
    xlabel('Fully-sampled Data')
    
    load(['accelfactor_',num2str(accelfactor),'/carotid_undersampled.mat']);
    mip = max(squeeze(rd_full(:,:,:,frames(ii))), [], 3);
    subplot(232)
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(frames(ii),2)),' ms'])
    xlabel('Under-sampled Data')
    ylabel(['acceleration factor = ',num2str(accelfactor)]);
    
    
    load(['accelfactor_',num2str(accelfactor),'/carotid_l1spirit.mat']);
    mip = max(squeeze(rd_full(:,:,:,frames(ii))), [], 3);
    subplot(233)
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(frames(ii),2)),' ms'])
    xlabel('L_1 SPIRiT')
    ylabel(['acceleration factor = ',num2str(accelfactor)]);
    
    
    load(['accelfactor_',num2str(accelfactor),'/carotid_distributed_cs_recon.mat']);
    mip = max(squeeze(rd_full(:,:,:,frames(ii))), [], 3);
    subplot(234)
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(frames(ii),2)),' ms'])
    xlabel('Distributed CS')
    ylabel(['acceleration factor = ',num2str(accelfactor)]);
    
    
    load(['accelfactor_',num2str(accelfactor),'/carotid_distributed_cs_recon_nocompensation.mat']);
    mip = max(squeeze(rd_full(:,:,:,frames(ii))), [], 3);
    subplot(235)
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(frames(ii),2)),' ms'])
    xlabel('Dist CS w/o Density Compensation')
    ylabel(['acceleration factor = ',num2str(accelfactor)]);
    
    
    load(['accelfactor_',num2str(accelfactor),'/carotid_distributed_cs_recon_phasecorrection.mat']);
    mip = max(squeeze(rd_full(:,:,:,frames(ii))), [], 3);
    subplot(236)
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(frames(ii),2)),' ms'])
    xlabel('Dist CS + Phase Correction')
    ylabel(['acceleration factor = ',num2str(accelfactor)]);
    
    
    pause(1)
    
    movie1(:,ii) = getframe(gcf); 
end

movie2avi(movie1,['accelfactor_',num2str(accelfactor),'/movies/carotid_dataset_complete.avi'],...
    'compression','None','fps',.75);