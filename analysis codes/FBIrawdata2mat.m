addpath(genpath('/home/mri-vh/Davi/functions'))

clear all;
close all;
clc

datasets = 5;
% datasets = [4 6 8:30];
dts = length(datasets);

% datasets [4 6 8:30]
% N = [112 256 10 16]; % [pe ro se coils]
% npefull = 160;

% dataset 5
N = [152 256 1 16]; % [pe ro se coils]
npefull = 160;

% N = [131 128 20 16]; % [pe ro se coils]
% npefull = N(1);
NAcq = 1;

rdc1 = zeros(npefull,N(2),N(4));
rdc3 = zeros(npefull,N(2),dts);

% figure; clf;
figr = floor(sqrt(dts));
figc = ceil(dts/figr);

for jj = 1
%     figure;
    for kk = 1:dts
        disp(['reading rawdata from dataset number ',num2str(datasets(kk))]),
        
            if datasets(kk) < 10
                raws2 = rawdataRead(['original raw data/Run1207.6904.0',num2str(datasets(kk)),'.raw.bin'], N(3), N(1), npefull, N(2), NAcq, N(4));
            else
                raws2 = rawdataRead(['original raw data/Run1207.6904.',num2str(datasets(kk)),'.raw.bin'], N(3), N(1), npefull, N(2), NAcq, N(4));
            end
        
        save(['Run1207.6904.',num2str(datasets(kk)),'_2.mat'],'raws2');
        
        if N(3) > 1
            
            for i = 1:N(4)
                rd1 = raws2(jj,:,:,i);
                rd1 = squeeze(rd1);
                rdc1(:,:,i) = rd1;
                
            end
            
            rdc2 = (abs(fftshift(fft2(fftshift(rdc1)))).^2);
            rdc3(:,:,kk) = sqrt(sum(rdc2,3));
            
            subplot(figr,figc,kk);
            imagesc(rot90(abs(rdc3(:,:,kk)),2)), title(['Image corresponding to dataset number n = ',num2str(datasets(kk))])
            colormap(gray)
        end
    end
end
