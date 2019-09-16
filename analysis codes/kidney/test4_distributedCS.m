addpath(genpath('/usr/local/MATLAB/R2013a/toolbox/sparseMRI_v0.2/'))
addpath(genpath('/home/mri-vh/Davi/functions'))

clear all;
% close all;
clc

% datasets = 9;
datasets = 6:13;
dts = length(datasets);

ncoils = 32;
N = [168, 128];
nslices = 24;

% pseudorandom sampling
% pdf = genPDF(N, 15, 5/8, 2, 0, 0);
% mask = genSampling(pdf,500,2);

load('mask.mat');

% CS parameters
TVWeight = 0.005;
xfmWeight = 0.01;
Itnlim = 15;

% Fourier sampling operator
FT = p2DFT(mask, N, 1, 2);
% XFM = Wavelet('Daubechies',8,4);
XFM = 1;

% CS initialization
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;

% Phase Correction
%phmask = zpad(hamming(6)*hamming(6)',N(1),N(2)); 	% mask to grab center frequency
%phmask = phmask/max(phmask(:));			 			% for low-order phase estimation and correction

figure; clf;
figr = floor(sqrt(dts));
figc = ceil(dts/figr);

rdc1 = zeros(N(1),N(2),ncoils);
rdc3 = zeros(N(1),N(2),dts);

im_dc = zeros(N(1),N(2),ncoils);
im_dc3 = zeros(N(1),N(2),dts);


for kk = 1:dts
    disp(['reconstructing dataset number ',num2str(datasets(kk))]),
    
    load(['Run2687.6904.',num2str(datasets(kk)),'.mat']);
    
    for slice = 1:nslices;
		%data = combinechannels(squeeze(raws2(slice,:,:,:)),3,2);
        %ph = exp(1i*angle((ifft2c(data.*phmask)))); % estimate phase for phase correction
        %FT = p2DFT(mask, N, ph, 2);
        %param.FT = FT;
        
        for i = 1:ncoils
            rd1a = raws2(slice,:,:,i);
            rd1a = squeeze(rd1a);
            rd2a = flipud(rd1a(1:floor((size(rd1a,1) + 1)/2),:));
            rd3a = (rd1a(floor((size(rd1a,1) + 1)/2 + 1):size(rd1a,1),:));
            rda = [rd2a; rd3a];
            rda = rda.*mask;
            
            clear rd1a rd2a rd3a
            
            % Methods:
            % 1st alternative: apply the Fourier transform, then compute the
            % difference, normalize the result using the pdf and perform the CS
            % recon;
            %         im_dc1 = FT'*rda - FT'*rdb;
            %         im_dc(:,:,i) = (im_dc1./pdf);
            %         rd = FT*im_dc1;
            %         rdc1(:,:,i) = rd;
            
            
            % 2nd alternative: apply the Fourier transform to the difference,
            % then use the pdf to normalize the data and finally perform the
            % recon.
            %         im_dc1 = FT'*(rda - rdb);
            %         im_dc(:,:,i) = (im_dc1./pdf);
            %         rd = FT*im_dc1;
            %         rdc1(:,:,i) = rd;
            
            % 3rd alternative (the one that presents better results until now):
            % compute the difference and normalize the signal using the pdf,
            % then apply the FT operator and perform the recon. This process
            % produced better angiograms (sparser) than the previous ones.
            rdiff = (rda);%./pdf;
            im_dc(:,:,i) = FT'*rdiff;
            rdc1(:,:,i) = rdiff;
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
        rdc3(:,:,slice) = sqrt(sum(rdc2,3));
        
        subplot(figr,figc,kk);
%         figure;
        imagesc(abs(rdc3(:,:,slice))), title(['slice = ', num2str(slice),', dataset number n = ',num2str(datasets(kk))])
        colormap(gray)
        
    end
    disp('Saving reconstructed dataset');
%    save(['kidney_dataset_',num2str(datasets(kk)),'_distributed_cs_recon.mat'], 'rdc3', 'mask','N','nslices');
	save(['kidney_dataset_',num2str(datasets(kk)),'_distributed_cs_no_compensated_recon.mat'], 'rdc3', 'mask','N','nslices');
end


