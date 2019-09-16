% addpath('/usr/MATLAB/R2013a/toolbox/sparseMRI_v0.2/')
% addpath('/usr/local/MATLAB/R2013a/toolbox/sparseMRI_v0.2/utils/')

clear all;
% close all;
clc

datasets = 12;
% datasets = 10:18;
dts = length(datasets);

ncoils = 15;
N = [128, 256];

% pseudorandom sampling
samplingfactor = 5/8;
pdf = genPDF(N, 20, samplingfactor, 2, 0, 1);
mask = genSampling(pdf,500,2);

% CS parameters
TVWeight = 0.05;
xfmWeight = 0.01;
Itnlim = 15;

% Fourier sampling operator
FT = p2DFT(mask, N, 1, 2);
XFM = Wavelet('Daubechies',4,4);
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


rdc1 = zeros(N(1),N(2),ncoils);
rdc3 = zeros(N(1),N(2),dts);


for kk = 1:dts
    disp(['reconstructing dataset number ',num2str(datasets(kk))]),
    
    load(['Run2730.6904.',num2str(datasets(kk)),'_2.mat']);
    
    for i = 1:ncoils
        rd1a = raws2(17,:,:,i);
        rd1a = squeeze(rd1a);
        
        rd2a = flipud(rd1a(1:floor((size(rd1a,1) + 1)/2),:));
        
        rd3a = (rd1a(floor((size(rd1a,1) + 1)/2 + 1):size(rd1a,1),:));
        rda = [rd2a; rd3a];
        rda = rda.*mask;
        
        rd1b = raws2(35,:,:,i);
        rd1b = squeeze(rd1b);
        
        rd2b = flipud(rd1b(1:floor((size(rd1b,1) + 1)/2),:));
        
        rd3b = (rd1b(floor((size(rd1b,1) + 1)/2 + 1):size(rd1b,1),:));
        rdb = [rd2b; rd3b];
        rdb = rdb.*mask;
        
        clear rd1a rd1b rd2a rd2b rd3a rd3b
        
        im_dc1 = FT'*rda - FT'*rdb;
        rd = FT*im_dc1;
        
        param.data = rd;
        
        im_dc = FT'*(rd./pdf);
        res = XFM*im_dc;
        
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

