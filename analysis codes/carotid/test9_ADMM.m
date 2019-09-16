addpath(genpath('/home/mri-vh/Davi/functions'))
addpath(genpath('/home/mri-vh/Davi/ADMM_CS_SENSE'))
addpath(genpath('/home/mri-vh/Davi/Marc_SPIRiT'))

clear all;
close all;
clc;


%% Data information
% datasets = 10;
% datasets = 11:18;
datasets = 10:18;
dts = length(datasets);

accelfactor = 10;
% mkdir(['accelfactor_',num2str(accelfactor)]);

load('params.mat');
params = [params 15];


%% Pseudorandom sampling scheme
load(['samplingmasks/sampling_mask_',num2str(accelfactor),'.mat']);
masku = repmat(mask,[1,1,params(4)]);


N = params(1);
N2 = params(1)*params(2);
NL = params(1)*params(2)*params(4);
% NL = params(2)*params(1);


% Total variance in two direction, a function version is also available
T1 = sparse([1:NL,2:NL],[1:NL,1:NL-1],[ones(1,NL),-ones(1,NL-1)]); % construct T1

T2 = sparse([1:NL,N2 + 1:NL],[1:NL,1:N2*(params(4)-1)],...
    [-ones(1,NL),ones(1,N2*(params(4)-1))]); % construct T2

% T2 = sparse([1:NL,params(2) + 1:NL],[1:NL,1:params(2)*(params(1)-1)],...
%     [-ones(1,NL),ones(1,params(2)*(params(1)-1))]); % construct T2

D = [T1;T2]; % construct tall matrix D

gama = 0.001; % dummy penalty, needs improvement
beta = 0.005;
opt.kx = params(1);
opt.ky = params(2);
opt.kz = params(3);
opt.coils = params(4);
opt.itermax = 10;
opt.tol = 1e-3;
opt.gama = 0.01;

wavWeight = 0.025;
nIter = 15;

kxkykzc = zeros(params(1),params(2),params(4),params(3));
rdc1 = zeros(params(1),params(2),opt.coils);
recon_data = zeros(params(1),params(2),params(4),params(3),dts);

for kk = 1:dts
    disp(['reconstructing dataset number ',num2str(datasets(kk))]),
    
    load(['Run2730.6904.',num2str(datasets(kk)),'_2.mat']);
    
    for slice = 1:params(3);
        disp(['working with slice number ',num2str(slice)])
        
        for i = 1:opt.coils
            rd1a = raws2(slice,:,:,i);
            rd1a = squeeze(rd1a);
            rd2a = flipud(rd1a(1:floor((size(rd1a,1) + 1)/2),:));
            rd3a = (rd1a(floor((size(rd1a,1) + 1)/2 + 1):size(rd1a,1),:));
            rda = [rd2a; rd3a];
            rda = rda.*mask;
            
            rd1b = raws2(slice + 18,:,:,i);
            rd1b = squeeze(rd1b);
            rd2b = flipud(rd1b(1:floor((size(rd1b,1) + 1)/2),:));
            rd3b = (rd1b(floor((size(rd1b,1) + 1)/2 + 1):size(rd1b,1),:));
            rdb = [rd2b; rd3b];
            rdb = rdb.*mask;
            
            clear rd1a rd1b rd2a rd2b rd3a rd3b
            
            rdiff = (rda - rdb)./pdf;
%             im_dc(:,:,i) = FT'*rdiff;
            rdc1(:,:,i) = rdiff;
        end
        
        b = gama*(ifft2c(masku.*rdc1)); % gama*DFT'*U'*kU to prepare for A'*b
        b = b(:);
        
        % set initial values
        xk = b;
        xk_p = xk;
        xk1 = zeros(NL,1);
        vk = zeros(2*NL,1); % dummy variable vk,uk
        iter = 1;
        itermax = 20;
        tol = 1e-4;
        
%         tic,
        while iter<itermax && norm(xk1-xk_p)>tol
            dk1 = shrink1(vk + D*xk, beta*gama);
            uk1 = vk + D*xk - dk1;
            [xk1,count] = CS_SENSE_recon(ones(1,1,opt.coils), masku, D, b + D'*(dk1 - uk1), xk, opt);
            
            % this version don't use the FT class
            
            xk_p = xk;
            xk = xk1;
            dk = dk1;
            vk = uk1;
            
%             imagesc(combinechannels(reshape(abs(xk),opt.kx,opt.ky,opt.coils),3,2)); colormap('gray');xlabel(num2str(slice)); drawnow
            iter = iter + 1;
        end
                
        kxkykzc(:,:,:,slice) = reshape(xk1,opt.kx,opt.ky,opt.coils);
        
%         slice
        
        ssx = 2^ceil(log2(opt.kx));
        ssy = 2^ceil(log2(opt.ky));
        ss = max(ssx, ssy);
        W = Wavelet('Daubechies',4,2);
        %W = Wavelet('Haar',2,3);
                
        x = fft2c(squeeze(kxkykzc(:,:,:,slice)));
        data = x;
        
        for n=1:nIter
            X = zpad(ifft2c(x),ss,ss,opt.coils);
            X = SoftThresh(W*(X),wavWeight);
            x = fft2c(crop(W'*(X),opt.kx,opt.ky,opt.coils)).*(1 - masku) + data;
        end
        
        iter = 1;
%         toc,
        kxkykzc(:,:,:,slice) = (x);
        
%         imgR = combinechannels(squeeze(ifft2c(x)),3,2); % reconstructed image
%         imagesc(abs(imgR)); colormap('gray');xlabel(num2str(slice)); drawnow
    end
    recon_data(:,:,:,:,kk) = kxkykzc;
end

% save(['ADMM_recon_data_us_',num2str(accelfactor),'_.mat'],'recon_data');
save(['ADMM_Wavelets_recon_data_us_',num2str(accelfactor),'_.mat'],'recon_data');

