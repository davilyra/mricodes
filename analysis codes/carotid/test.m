clear all;
% close all;
clc

reducedfov = 1; % 1 = use reduced fov; 0 = use regular fov;

datasets = 13;
% datasets = 10:18;
dts = length(datasets);

ncoils = 15;

figure; clf;
figr = floor(sqrt(dts));
figc = ceil(dts/figr);

if reducedfov == 1
    rdc1 = zeros(128,256,ncoils);
    rdc3 = zeros(128,256,dts);
else
    rdd = zeros(256,256,ncoils);
    rdc3 = zeros(256,256,dts);
end

for jj = 1
    figure;
    for kk = 1:dts
        disp(['reconstructing dataset number ',num2str(datasets(kk))]),
        
%             if datasets(kk) < 10
%                 raws2 = rawdataRead(['Run2730.6904.0',num2str(datasets(kk)),'.raw.bin'], 36, 128, 128, 256, 1, ncoils);
%             else
%                 raws2 = rawdataRead(['Run2730.6904.',num2str(datasets(kk)),'.raw.bin'], 36, 128, 128, 256, 1, ncoils);
%             end
        
%         load(['Run2730.6904.',num2str(datasets(kk)),'.mat']);
        save(['Run2730.6904.',num2str(datasets(kk)),'_2.mat'],'raws2');
        
        for i = 1:ncoils
            rd1 = raws2(jj,:,:,i);
            rd1 = squeeze(rd1);
            
            rd2 = flipud(rd1(1:floor((size(rd1,1) + 1)/2),:));
            
            rd3 = (rd1(floor((size(rd1,1) + 1)/2 + 1):size(rd1,1),:));
            rd = [rd2; rd3];
            
            if reducedfov == 1
                rdc1(:,:,i) = rd;
                %             figure, imagesc((abs(squeeze(rdc1(:,:,i))))),colormap(gray)
                %             figure, imagesc((abs(squeeze(fftshift(fft2(fftshift(rdc1(:,:,i)))))))),colormap(gray)
            else
                rdd(1:2:end, :, i) = rd;
                %             figure, imagesc((abs(squeeze(rdd(:,:,i))))),colormap(gray)
                %             figure, imagesc((abs(squeeze(fftshift(fft2(fftshift(rdd(:,:,i)))))))),colormap(gray)
            end
            
        end
        
        if reducedfov == 1
            %     Using reduced field-of-view reconstruction
            rdc2 = (abs(fftshift(fft2(fftshift(rdc1)))).^2);
            rdc3(:,:,kk) = sqrt(sum(rdc2,3));
            
            subplot(figr,figc,kk);
            imagesc(rot90(abs(rdc3(:,:,kk)),2)), title(['Image corresponding to dataset number n = ',num2str(datasets(kk))])
            colormap(gray)
        else
            %     Reconstruction with the entire field-of-view
            rdc2 = (abs(fftshift(fft2(fftshift(rdd)))).^2);
            rdc3(:,:,kk) = sqrt(sum(rdc2,3));
            
            subplot(figr,figc,kk);
            imagesc(rot90(abs(rdc3(:,:,kk)),2)), title(['Image corresponding to dataset number n = ',num2str(datasets(kk))])
            colormap(gray)
        end
    end
end
