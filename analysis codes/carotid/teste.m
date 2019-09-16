clear all;
% close all;
clc

% rdd = zeros(188,128,32);
% rdd = zeros(336,128,32);

for kk = 13:2:27
    disp(['reconstructing data from d2362368/128/769ataset number ',num2str(kk)]),
    for i = 1:32
        %     raws2(:,:,:,i) = rawdataRead(['Run2687.6904.5.',num2str(i)], 24, 92, 92, 128, 1, 1, 1);
        % for the data contained in .5, the dimensions [24, 94, 128, 32] are the correct
        
        raws2(:,:,:,i) = rawdataRead(['Run2687.6904.',num2str(kk),'.',num2str(i)], 24, 168, 168, 128, 1, 1, 1);
        
        rd1 = raws2(16,:,:,i);
        rd1 = squeeze(rd1);
        
        rd2 = flipud(rd1(1:floor((end + 1)/2),:));
        %     rd2 = (rd1(1:floor((end +1)/2),:));
        
        rd3 = (rd1(floor((end + 1)/2 + 1):end,:));
        rd = [rd2; rd3];
        
        % 	rdd(1:2:end, :, i) = rd;
        rdc1(:,:,i) = rd;
        
        % 	figure, imagesc((abs(squeeze(rdc1(:,:,i))))),colormap(gray)
        % 	figure, imagesc((abs(squeeze(fftshift(fft2(fftshift(rdc1(:,:,i)))))))),colormap(gray)
        
        
        % 	figure, imagesc((abs(squeeze(rdd(:,:,i))))),colormap(gray)
        % 	figure, imagesc((abs(squeeze(fftshift(fft2(fftshift(rdd(:,:,i)))))))),colormap(gray)
        
    end
    
    rdc2 = (abs(fftshift(fft2(fftshift(rdc1)))).^2);
    % rdc2 = (abs(fftshift(fft2(fftshift(rdd)))).^2);
    
    rdc3 = sqrt(sum(rdc2,3));
    
    figure,imagesc(abs(rdc3)), title(['Image corresponding to dataset number n = ',num2str(kk)])
end