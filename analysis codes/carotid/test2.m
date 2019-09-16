clear all;
% close all;
clc

reducedfov = 1; % 1 = use reduced fov; 0 = use regular fov;

datasets = 13;
% datasets = 10:18;
dts = length(datasets);

ncoils = 15;
N = [128, 256];

% pseudorandom sampling
pdf = genPDF(N, 15, 5/8, 2, 0, 0);
mask = genSampling(pdf,500,2);

% CS parameters
TVWeight = 0.05;
xfmWeight = 0.01;
Itnlim = 25;

% Fourier sampling operator
FT = p2DFT(mask, N, 1, 2);
XFM = Wavelet('Daubechies',8,4);
% XFM = 1;

param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;


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

for kk = 1:dts
    disp(['reconstructing dataset number ',num2str(datasets(kk))]),
    
    load(['Run2730.6904.',num2str(datasets(kk)),'.mat']);
    
    for i = 1:ncoils
        rd1 = raws2(17,:,:,i);
        rd1 = squeeze(rd1);
        
        rd2 = flipud(rd1(1:floor((size(rd1,1) + 1)/2),:));
        
        rd3 = (rd1(floor((size(rd1,1) + 1)/2 + 1):size(rd1,1),:));
        rd = [rd2; rd3];
        rd = rd.*mask;
        
        param.data = rd;
        
        im_dc = FT'*(rd./pdf);
        res = XFM*im_dc;
        
        for n=1:5
            res = fnlCg(res,param);
%             im_res = XFM'*res;
%             figure(100), imshow(abs(im_res),[]), drawnow
        end
        
        if reducedfov == 1
            rdc1(:,:,i) = XFM'*res;
%             figure, imagesc((abs(squeeze(rdc1(:,:,i))))),colormap(gray)
%             figure, imagesc((abs(squeeze(fftshift(fft2(fftshift(rdc1(:,:,i)))))))),colormap(gray)
        else
            rdd(1:2:end, :, i) = XFM'*res;
%             figure, imagesc((abs(squeeze(rdd(:,:,i))))),colormap(gray)
%             figure, imagesc((abs(squeeze(fftshift(fft2(fftshift(rdd(:,:,i)))))))),colormap(gray)
        end
        
    end
    
    if reducedfov == 1
        %     Using reduced field-of-view reconstruction
%         rdc2 = (abs(fftshift(fft2(fftshift(rdc1)))).^2);
        rdc2 = (abs(rdc1).^2);
        rdc3(:,:,kk) = sqrt(sum(rdc2,3));
        
        subplot(figr,figc,kk);
        imagesc(abs(rdc3(:,:,kk))), title(['Image corresponding to dataset number n = ',num2str(datasets(kk))])
        colormap(gray)
    else
        %     Reconstruction with the entire field-of-view
        rdc2 = (abs(rdd).^2);
        rdc3(:,:,kk) = sqrt(sum(rdc2,3));
        
        subplot(figr,figc,kk);
        imagesc(abs(rdc3(:,:,kk))), title(['Image corresponding to dataset number n = ',num2str(datasets(kk))])
        colormap(gray)
    end
end

