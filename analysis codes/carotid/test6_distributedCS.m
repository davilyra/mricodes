clear all;
% close all;
clc

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
Itnlim = 15;

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

rdc1 = zeros(128,256,ncoils);
im_dc = zeros(128,256,ncoils);
rdc3 = zeros(128,256,dts);
im_dc3 = zeros(128,256,dts);

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
        
        rdc1(:,:,i) = rd;
        im_dc(:,:,i) = FT'*(rdc1(:,:,i)./pdf);
    end
    
    im_dc2 = (abs(im_dc).^2);
    im_dc3(:,:,kk) = sqrt(sum(im_dc2,3));
    
	res = XFM*im_dc3(:,:,kk);

    
    for i = 1:ncoils
        
        param.data = rdc1(:,:,i);
        
        for n=1:5
            res = fnlCg(res,param);
        end
        
        rdc1(:,:,i) = XFM'*res;
    end
    
    
    rdc2 = (abs(rdc1).^2);
    rdc3(:,:,kk) = sqrt(sum(rdc2,3));
    
    subplot(figr,figc,kk);
    imagesc(abs(rdc3(:,:,kk))), title(['Image corresponding to dataset number n = ',num2str(datasets(kk))])
    colormap(gray)
    
end

