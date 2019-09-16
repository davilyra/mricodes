%% Organizing workspace
clear all; close all; clc;

% Sorting the data according to the BBTI
folder2 = [400 500 600 700 800 900 300 1000 1100]; % BBTI values from the performed study
folder1 = [10 11 12 13 14 15 16 17 18]; % Datasets id
folder = [folder1' folder2'];
folder = sortrows(folder,2);
clear folder1 folder2

accelfactor = 8;

load('/home/mri-vh/Davi/carotid/params.mat')
mkdir(['accelfactor_',num2str(accelfactor),'/movies']);


rd_full = zeros(params(2),params(1),params(3),size(folder,1));
movie1 = moviein(10);

for i = 1:size(folder,1)
    load(['accelfactor_',num2str(accelfactor),'/carotid_dataset_',num2str(folder(i)),'_fullysampled.mat']);

    for jj = 1:params(3)
        rd_full(:,:,jj,i) = rot90(rdc3(:,:,jj),-1);
    end
    
    mip = max(squeeze(rd_full(:,:,:,i)), [], 3);
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(i,2)),' ms'])
    xlabel('Fully-sampled Data')
    pause(1)
    
    movie1(:,i) = getframe(gcf);
end
save(['accelfactor_',num2str(accelfactor),'/carotid_fullysampled.mat'],'rd_full');
movie2avi(movie1,['accelfactor_',num2str(accelfactor),'/movies/carotid_fullysampled.avi'],...
    'compression','None','fps',.75);


for i = 1:size(folder,1)
    load(['accelfactor_',num2str(accelfactor),'/carotid_dataset_',num2str(folder(i)),'_undersampled.mat']);

    for jj = 1:params(3)
        rd_full(:,:,jj,i) = rot90(rdc3_us(:,:,jj),-1);
    end
    
    mip = max(squeeze(rd_full(:,:,:,i)), [], 3);
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(i,2)),' ms'])
    xlabel('Undersampled Data')
    pause(1)
    
    movie1(:,i) = getframe(gcf);
end
save(['accelfactor_',num2str(accelfactor),'/carotid_undersampled.mat'],'rd_full');
movie2avi(movie1,['accelfactor_',num2str(accelfactor),'/movies/carotid_undersampled.avi'],...
    'compression','None','fps',.75);


for i = 1:size(folder,1)
    load(['accelfactor_',num2str(accelfactor),'/carotid_dataset_',num2str(folder(i)),'_distributed_cs_recon.mat']);
  
    for jj = 1:params(3)
        rd_full(:,:,jj,i) = rot90(rdc3(:,:,jj),1);
    end
    
    mip = max(squeeze(rd_full(:,:,:,i)), [], 3);
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(i,2)),' ms'])
    xlabel('Distributed CS')
    pause(1)
    
    movie1(:,i) = getframe(gcf);
end
save(['accelfactor_',num2str(accelfactor),'/carotid_distributed_cs_recon.mat'],'rd_full');
movie2avi(movie1,['accelfactor_',num2str(accelfactor),'/movies/carotid_distributed_cs_recon.avi'],...
    'compression','None','fps',.75);


for i = 1:size(folder,1)
    load(['accelfactor_',num2str(accelfactor),'/carotid_dataset_',num2str(folder(i)),'_distributed_cs_recon_nocompensation.mat']);
    for jj = 1:params(3)
        rd_full(:,:,jj,i) = rot90(rdc3(:,:,jj),1);
    end
    
	mip = max(squeeze(rd_full(:,:,:,i)), [], 3);
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(i,2)),' ms'])
    xlabel('Dist CS w/o Density Compensation')
    pause(1)
    
    movie1(:,i) = getframe(gcf);
end
save(['accelfactor_',num2str(accelfactor),'/carotid_distributed_cs_recon_nocompensation.mat'],'rd_full');
movie2avi(movie1,['accelfactor_',num2str(accelfactor),'/movies/carotid_distributed_cs_recon_nocompensation.avi'],...
    'compression','None','fps',.75);


for i = 1:size(folder,1)
    load(['accelfactor_',num2str(accelfactor),'/carotid_dataset_',num2str(folder(i)),'_distributed_cs_recon_phasecorrection.mat']);

    for jj = 1:params(3)
        rd_full(:,:,jj,i) = rot90(rdc3(:,:,jj),1);
    end
    
    mip = max(squeeze(rd_full(:,:,:,i)), [], 3);
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(i,2)),' ms'])
    xlabel('Dist CS + Phase Correction')
    pause(1)
    
    movie1(:,i) = getframe(gcf);
end
save(['accelfactor_',num2str(accelfactor),'/carotid_distributed_cs_recon_phasecorrection.mat'],'rd_full');
movie2avi(movie1,['accelfactor_',num2str(accelfactor),'/movies/carotid_distributed_cs_recon_phasecorrection.avi'],...
    'compression','None','fps',.75);


for i = 1:size(folder,1)
    load(['accelfactor_',num2str(accelfactor),'/carotid_dataset_',num2str(folder(i)),'_l1spirit_recon.mat']);

    for jj = 1:params(3)
        rd_full(:,:,jj,i) = rot90(rdc3(:,:,jj),1);
    end
    
    mip = max(squeeze(rd_full(:,:,:,i)), [], 3);
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(i,2)),' ms'])
    xlabel('L_1 SPIRiT')
    pause(1)
    
    movie1(:,i) = getframe(gcf);
end
save(['accelfactor_',num2str(accelfactor),'/carotid_l1spirit.mat'],'rd_full');
movie2avi(movie1,['accelfactor_',num2str(accelfactor),'/movies/carotid_l1spirit_recon.avi'],...
    'compression','None','fps',.75);

load(['ADMM_recon_data_us_',num2str(accelfactor),'.mat']);
for i = 1:size(folder,1)    
    imgR = combinechannels(squeeze(recon_data(:,:,:,:,i)),3,2);
    mip = max(squeeze(imgR(:,:,:,:)), [], 3);
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(i,2)),' ms'])
    xlabel('ADMM CS')
    pause(1)
    
    movie1(:,i) = getframe(gcf);
end
save(['accelfactor_',num2str(accelfactor),'/carotid_l1spirit.mat'],'rd_full');
movie2avi(movie1,['accelfactor_',num2str(accelfactor),'/movies/carotid_admm_cs.avi'],...
    'compression','None','fps',.75);

load(['IRLS_recon_data_us_',num2str(accelfactor),'.mat']);
for i = 1:size(folder,1)    
    imgR = combinechannels(squeeze(recon_data(:,:,:,:,i)),3,2);
    mip = max(squeeze(imgR(:,:,:,:)), [], 3);
    imshow(abs(mip),[]),colormap(gray)
    title(['BBTI = ',num2str(folder(i,2)),' ms'])
    xlabel('IRLS CS')
    pause(1)
    
    movie1(:,i) = getframe(gcf);
end
save(['accelfactor_',num2str(accelfactor),'/carotid_l1spirit.mat'],'rd_full');
movie2avi(movie1,['accelfactor_',num2str(accelfactor),'/movies/carotid_irls_cs.avi'],...
    'compression','None','fps',.75);