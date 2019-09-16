%% undersampling the data

% preparing the workspace
clear all; close all; clc;

load('params.mat');

accelfactor = 2;

app_select = 2; %[1] = full; [2] = fbi;
app_select2 = 1; %[1] = pese; [2] = pero;


if app_select2 == 1
    application = 'FBI/pese';
    
    mkdir([application(5:end),'/accelfactor_',num2str(accelfactor)]);
    load(['/home/mri-vh/Davi/samplingmasks/FBI/pese/sampling_mask_',num2str(accelfactor),'.mat']);
    mask = permute(repmat(mask,[1,1,params(2),params(4),params(5)]),[3 1 2 4 5]);
    pdf = permute(repmat(pdf,[1,1,params(2),params(4),params(5)]),[3 1 2 4 5]);
    
elseif app_select2 == 2
    application = 'FBI/pero';
    
    mkdir([application(5:end),'/accelfactor_',num2str(accelfactor)]);
    load(['/home/mri-vh/Davi/samplingmasks/FBI/pero/sampling_mask_',num2str(accelfactor),'.mat']);
    mask = permute(repmat(mask,[1,1,params(3),params(4),params(5)]),[2 1 3 4 5]);
    pdf = permute(repmat(pdf,[1,1,params(3),params(4),params(5)]),[2 1 3 4 5]);
    
end

% save(['/home/mri-vh/Davi/samplingmasks/',application,'/cs_sampling_mask_',num2str(accelfactor),'.mat'],'mask','pdf','-v7.3');


switch app_select
    case 1
        
        load('rawdata_full.mat');
        rd_us = rd_full.*mask./pdf;
        
        save([application(5:end),'/accelfactor_',num2str(accelfactor),...
            '/rawdata_',application(5:end),'_us_',num2str(accelfactor),'.mat'],'rd_us','-v7.3');
        
    case 2
        
        load('rawdata_fbi.mat');
        rd_fbi_us = rd_fbi.*mask;
        
        save([application(5:end),'/accelfactor_',num2str(accelfactor),...
            '/rawdata_fbi1_',application(5:end),'_us_',num2str(accelfactor),'.mat'],'rd_fbi_us','-v7.3');
        
%         rd_fbi_us = rd_fbi.*mask./pdf;
%         
%         save([application(5:end),'/accelfactor_',num2str(accelfactor),...
%             '/rawdata_fbi_',application(5:end),'_us_',num2str(accelfactor),'.mat'],'rd_fbi_us','-v7.3');
        
end