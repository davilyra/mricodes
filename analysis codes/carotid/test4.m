clear all;
% close all;
clc

% datasets = 10;
datasets = 10:18;
dts = length(datasets);

ncoils = 15;
N = [128, 256];
nslices = 18;
accelfactor = 8;
mkdir(['accelfactor_',num2str(accelfactor)]);

% figure; clf;
figr = floor(sqrt(dts));
figc = ceil(dts/figr);


% rdc1_a = zeros(128,256,ncoils);
% rdc1_b = zeros(128,256,ncoils);

rdiff = zeros(128,256,ncoils);
rdiff_us = zeros(128,256,ncoils);
rdc3 = zeros(128,256,nslices);
rdc3_us = zeros(128,256,nslices);

% rdc3_a = zeros(128,256,dts);
% rdc3_b = zeros(128,256,dts);
load(['samplingmasks/sampling_mask_',num2str(accelfactor),'.mat']);

for kk = 1:dts
    disp(['reconstructing dataset number ',num2str(datasets(kk))]),
    
    load(['Run2730.6904.',num2str(datasets(kk)),'_2.mat']);

    for slice = 1:nslices;
        for i = 1:ncoils
            rd1 = raws2(slice,:,:,i);
            rd1 = squeeze(rd1);
            rd2 = flipud(rd1(1:floor((size(rd1,1) + 1)/2),:));
            rd3 = (rd1(floor((size(rd1,1) + 1)/2 + 1):size(rd1,1),:));
            rd_a = [rd2; rd3];
            rda_us = rd_a.*mask;
            
            rd1 = raws2(slice + 18,:,:,i);
            rd1 = squeeze(rd1);
            rd2 = flipud(rd1(1:floor((size(rd1,1) + 1)/2),:));
            rd3 = (rd1(floor((size(rd1,1) + 1)/2 + 1):size(rd1,1),:));
            rd_b = [rd2; rd3];
            rdb_us = rd_b.*mask;
            
            rdiff(:,:,i) = rd_a - rd_b;
            rdiff_us(:,:,i) = rda_us - rdb_us;
            
        end
        
        rdiff_c = (abs(fftshift(fft2(fftshift(rdiff)))).^2);
        rdc3(:,:,slice) = sqrt(sum(rdiff_c,3));
        
        rdiff_c_us = (abs(fftshift(fft2(fftshift(rdiff_us)))).^2);
        rdc3_us(:,:,slice) = sqrt(sum(rdiff_c_us,3));
        
%         rdc2_a = (abs(fftshift(fft2(fftshift(rdc1_a)))).^2);
%         rdc2_b = (abs(fftshift(fft2(fftshift(rdc1_b)))).^2);
%         rdc3_a(:,:,kk) = sqrt(sum(rdc2_a,3));
%         rdc3_b(:,:,kk) = sqrt(sum(rdc2_b,3));
%         
    end
    save(['accelfactor_',num2str(accelfactor),'/carotid_dataset_',num2str(datasets(kk)),'_fullysampled.mat'], 'rdc3', 'N','nslices');
    save(['accelfactor_',num2str(accelfactor),'/carotid_dataset_',num2str(datasets(kk)),'_undersampled.mat'], 'rdc3_us', 'N','nslices');

end
