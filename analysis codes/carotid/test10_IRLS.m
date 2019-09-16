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

accelfactor = 8;
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
T = diag(ones(N,1),0) + diag(-ones(N-1,1),-1);

A = sparse(1:N,1:N,ones(N,1)); % A is the identity matrix

itermax = 100; % iteration parameters
itermax2 = 100;
tol1 = 1e-4;
tol2 = 1e-3;

L = 100;
lambda = 0.01;
beta = 0.02;
mu = 1e-6;


wavWeight = 0.0025;
nIter = 15;

kxkykzc = zeros(params(1),params(2),params(4),params(3));
rdc1 = zeros(params(1),params(2),params(4));
xr = zeros(params(1),params(2));
recon_data = zeros(params(1),params(2),params(4),params(3),dts);

for kk = 1:dts
    disp(['reconstructing dataset number ',num2str(datasets(kk))]),
    
    load(['Run2730.6904.',num2str(datasets(kk)),'_2.mat']);
    
    for slice = 1:params(3);
        disp(['working with slice number ',num2str(slice)])
        
        for i = 1:params(4)
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
            rdc1(:,:,i) = rdiff;
            
            b = (ifft2c(mask.*rdc1(:,:,i))); % gama*DFT'*U'*kU to prepare for A'*b
            
%         tic           
            for ii = 1:params(2)
                b1 = b(:,ii); % treat each column individually
                
                x_k = randn(N,1); % set initial values, need to set up
                x_k1 = zeros(N,1);
                x_kp = x_k;
                iter1 = 1;
                
                while iter1<itermax && norm(x_k1 - x_kp)>tol1
                    zk = min(L,1./(2*abs(T*x_k))); % get weighted zk from each xk
                    Wzk = sparse(1:N,1:N,zk); % get diagnol matrix Wzk
                    AL = [A; sqrt(lambda)*Wzk*T]; % get the tall matrix AL
                    Q = AL'*AL; % get matrix Q in CG
                    [x_k1,count] = CG_Q(Q, b1, x_k, itermax2, tol2);
                    % implement CG method, the previous value x_k is used to converge faster
                    
                    x_kp = x_k; % update varibale
                    x_k = x_k1;
                    iter1 = iter1+1;
                end
                
                xr(:,ii) = x_k;
            end
%             imshow(abs(xr),[])
%             pause(.5)
            kxkykzc(:,:,i,slice) = fft2c(xr);
        end
        
        ssx = 2^ceil(log2(params(1)));
        ssy = 2^ceil(log2(params(2)));
        ss = max(ssx, ssy);
        W = Wavelet('Daubechies',4,4);
        
%         x = ifft2c(squeeze(kxkykzc(:,:,:,slice)));
        x = (squeeze(kxkykzc(:,:,:,slice)));
        data = x;
              
        for n=1:nIter
            X = zpad(ifft2c(x),ss,ss,params(4));
            X = SoftThresh(W*(X),wavWeight);
            x = fft2c(crop(W'*(X),params(1),params(2),params(4))).*(1 - masku) + data;
        end
       
%         toc,
        kxkykzc(:,:,:,slice) = (x);
        
%         imgR = combinechannels(squeeze(ifft2c(x)),3,2); % reconstructed image
%         imagesc(abs(imgR)); colormap('gray');xlabel(num2str(slice)); drawnow
    end
    recon_data(:,:,:,:,kk) = kxkykzc;
end

% save(['IRLS_recon_data_us_',num2str(accelfactor),'_.mat'],'recon_data');
save(['IRLS_Wavelets_recon_data_us_',num2str(accelfactor),'_.mat'],'recon_data');
