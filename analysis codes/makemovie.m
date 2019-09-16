addpath(genpath('/home/mri-vh/Davi/functions'))

%% Clearing the Workspace
clear all; close all; clc

application = 'pese';

accelfactor = 8; % CS reduction factor
frames = [1:25 25];
nframes = length(frames); % number of BBTI's

% folder2 = [400 500 600 700 800 900 300 1000 1100]; % BBTI values from the performed study
% folder1 = [10 11 12 13 14 15 16 17 18]; % Datasets id
% folder = [folder1' folder2'];
% folder = sortrows(folder,2);
% clear folder1 folder2

load('/home/mri-vh/Davi/FBI Data/1207 data/params.mat')

load('rawdata_fbi');
load([application,'/accelfactor_',num2str(accelfactor),'/rawdata_fbi1_',application,'_us_',num2str(accelfactor),'.mat']);
load([application,'/accelfactor_',num2str(accelfactor),'/admm_cs_recon.mat']); kxkykzc1 = kxkykzc;
load([application,'/accelfactor_',num2str(accelfactor),'/irls_cs_recon1.mat']);
% load([application,'/accelfactor_',num2str(accelfactor),'/fbi_l1_spirit_recon.mat']); recon_data = recon_data1;
load([application,'/accelfactor_',num2str(accelfactor),'/fbi_s-spirit_recon.mat']);

movie1 = moviein(nframes);

for ii = 1:nframes;
    mip = max(squeeze(combinechannels(fft2c(rd_fbi(:,:,:,:,frames(ii))),4,2)), [], 3);
    subplot(231)
%     subplot(221)
    imshow(abs(mip),[]),colormap(gray)
    title(['time delay number ',num2str(frames(ii))]);
    xlabel('Fully-sampled Data')
    
    mip = max(squeeze(combinechannels(fft2c(rd_fbi_us(:,:,:,:,frames(ii))),4,2)), [], 3);
    subplot(232)
%     subplot(222)
    imshow(abs(mip),[]),colormap(gray)
    xlabel('Zero-filled Data')
    ylabel(['acceleration factor = ',num2str(accelfactor)]);
    
   
    mip = max(squeeze(combinechannels(fft2c(kxkykzc1(:,:,:,:,frames(ii))),4,2)), [], 3);
    subplot(233)
%     subplot(223)
    imshow(fliplr(abs(mip)),[]),colormap(gray)
%     imshow(rot90(abs(mip),2),[]),colormap(gray)
    xlabel('ADMM Data')
    ylabel(['acceleration factor = ',num2str(accelfactor)]);
    
    mip = max(squeeze(combinechannels(fft2c(kxkykzc(:,:,:,:,frames(ii))),4,2)), [], 3);
    subplot(234)
%     subplot(224)
    imshow(fliplr(abs(mip)),[]),colormap(gray)
    xlabel('IRLS Data')
    ylabel(['acceleration factor = ',num2str(accelfactor)]);

    mip = max(squeeze(combinechannels(fft2c(recon_data1(:,:,:,:,frames(ii))),4,2)), [], 3);
    subplot(235)
%     subplot(224)
    imshow(abs(mip),[]),colormap(gray)
    xlabel('L1 s-SPIRiT Data')
    ylabel(['acceleration factor = ',num2str(accelfactor)]);
    
%     mip = max(squeeze(combinechannels(fft2c(recon_data(:,:,:,:,frames(ii))),4,2)), [], 3);
%     subplot(236)
% %     subplot(224)
%     imshow(abs(mip),[]),colormap(gray)
%     xlabel('L1 SPIRiT Data')
%     ylabel(['acceleration factor = ',num2str(accelfactor)]);
    
    
    pause(1)
    
    movie1(:,ii) = getframe(gcf); 
end

movie2avi(movie1,[application,'_accelfactor_',num2str(accelfactor),'_movie_complete.avi'],...
    'compression','None','fps',.75);