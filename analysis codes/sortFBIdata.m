% Prepare Workspace
clear; close all; clc;

% Data Information
BBTI = [370 350 390:20:830];
dataNum = [6 4 8:30];
BBTIdataNum = [dataNum' BBTI'];

params = [160 256 10 16 size(BBTIdataNum,1)]; % [pe ro se coils nBBTIs];

type = 2; %[1] = full data; [2] = fbi data;

switch type
    case 1
        rd_full = zeros(params(2),params(1),params(3), params(4), params(5));
        
        % Combining the data into one single tensor
        for kk = 1:size(BBTIdataNum,1)
            load(['Run1207.6904.',num2str(BBTIdataNum(kk)),'_2.mat']); %  size(raws2) = [nse npe nro ncoils]
            rd_full(:,:,:,:,kk) = permute(raws2,[3 2 1 4]); % size(rd_full) = [nro npe nse ncoils nBBTIs]
        end
        
        save('rawdata_full.mat','rd_full','-v7.3');
        
    case 2
        rd_fbi = zeros(params(2),params(1),params(3), params(4), params(5));
        
        load('Run1207.6904.6_2.mat');
        raws1 = permute(raws2,[3 2 1 4]);
        
        % Combining the data into one single tensor
        for kk = 1:size(BBTIdataNum,1)
            load(['Run1207.6904.',num2str(BBTIdataNum(kk)),'_2.mat']); %  size(raws2) = [nse npe nro ncoils]
            rd_fbi(:,:,:,:,kk) = permute(raws2,[3 2 1 4]); % size(rd_fbi) = [nro npe nse ncoils nBBTIs]
            rd_fbi(:,:,:,:,kk) = rd_fbi(:,:,:,:,kk) - raws1; % obtaining the angiograms by subtraction
        end
        
        save('rawdata_fbi.mat','rd_fbi','-v7.3');
        
end

% save('params.mat','params');
% save('BBTIsort.mat','BBTIdataNum');
