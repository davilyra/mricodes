%% Analysis of the reconstructed datasets
clear alll; close all; clc

accelfactor = 8;
cd(['accelfactor_',num2str(accelfactor)]);
mkdir('images')

datasets = 10:18;
dts = length(datasets);

for kk = 1:dts
    load(['carotid_dataset_',num2str(datasets(kk)),'_fullysampled.mat'])
    mip = rot90(max(rdc3, [], 3),-1);
    mip = mip/max(mip(:));
    imwrite(mip,['images/dataset_',num2str(datasets(kk)),'_mip_fs.png'],'PNG');
    
    load(['carotid_dataset_',num2str(datasets(kk)),'_distributed_cs_recon.mat'])
    mip = rot90(max(rdc3, [], 3),1);
    mip = mip/max(mip(:));
    imwrite(mip,['images/dataset_',num2str(datasets(kk)),'_mip_distcs.png'],'PNG');
    
    load(['carotid_dataset_',num2str(datasets(kk)),'_distributed_cs_recon_phasecorrection.mat'])
    mip = rot90(max(rdc3, [], 3),1);
    mip = mip/max(mip(:));
    imwrite(mip,['images/dataset_',num2str(datasets(kk)),'_mip_distcs_phasecorrection.png'],'PNG');
    
    load(['carotid_dataset_',num2str(datasets(kk)),'_distributed_cs_recon_nocompensation.mat'])
    mip = rot90(max(rdc3, [], 3),1);
    mip = mip/max(mip(:));
    imwrite(mip,['images/dataset_',num2str(datasets(kk)),'_mip_distcs_nocompensation.png'],'PNG');
    
    load(['carotid_dataset_',num2str(datasets(kk)),'_undersampled.mat'])
    mip = rot90(max(rdc3_us, [], 3),-1);
    mip = mip/max(mip(:));
    imwrite(mip,['images/dataset_',num2str(datasets(kk)),'_mip_us.png'],'PNG');
    
    load(['carotid_dataset_',num2str(datasets(kk)),'_l1spirit_recon.mat'])
    mip = rot90(max(rdc3, [], 3),1);
    mip = mip/max(mip(:));
    imwrite(mip,['images/dataset_',num2str(datasets(kk)),'_mip_l1spirit.png'],'PNG');
    
end