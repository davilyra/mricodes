% denoise based on single column
% use IRLS to solve l1 norm

addpath(genpath('/home/mri-vh/Davi/functions'))
addpath(genpath('/home/mri-vh/Davi/ADMM_CS_SENSE'))

clear all;
% close all;
clc;

% loading data size parameters
load('params.mat');

% recon parameters
N = params(1);
A = sparse(1:N,1:N,ones(N,1)); % A is the identity matrix

itermax = 100; % iteration parameters
itermax2 = 100;
tol1 = 1e-4;
tol2 = 1e-3;

L = 100;
lambda = 0.01;

T = diag(ones(N,1),0) + diag(-ones(N-1,1),-1); % total variance matrix in kx direction

xr = zeros(params(1),params(3));
kxkykzc = zeros(params(2),params(1),params(3),params(4),params(5));

application = 'FBI/pese';

% data and recon information
% accelfactors = [8 10];
accelfactors = [2 4 5 8 10];
uslth = length(accelfactors);

for iii = 1:uslth
    accelfactor = accelfactors(iii);
    
    
    mkdir([application(5:end),'/accelfactor_',num2str(accelfactor)]);
    
    % load(['/home/mri-vh/Davi/samplingmasks/',application,'/sampling_mask_',...
    %     num2str(accelfactor),'.mat']);
    
    load([application(5:end),'/accelfactor_',num2str(accelfactor),...
        '/rawdata_fbi1_',application(5:end),'_us_',num2str(accelfactor),'.mat']);
    
    for bbti = 2:params(5);
        % for bbti = 15;
        tic,
        for kx = 1:params(2)
            img = fftshift(ifft2(fftshift(squeeze(rd_fbi_us(kx,:,:,:,bbti)))));
            
            kx
            
            for coils = 1:params(4)
                b = squeeze(img(:,:,coils)); % the measured image as b
                
                for i = 1:params(3)
                    b1 = b(:,i); % treat each column individually
                    
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
                        %         implement CG method, the previous value x_k is used to converge faster
                        
                        x_kp = x_k; % update varibale
                        x_k = x_k1;
                        iter1 = iter1+1;
                    end
                    
                    xr(:,i) = x_k;
                    
                end
                
                kxkykzc(kx,:,:,coils,bbti) = ifft2c(xr);
                
            end
            
        end
    end
    
    % masku = repmat(mask,[1,1,params(4),params(2),params(5)]);
    % masku = squeeze(permute(masku(:,:,1,:,15),[4 1 3 2 5]));
    
    % nIter = 15;
    % wavWeight = 0.025;
    
    % x = ifft2c((kxkykzc));
    
    % data = x;
    %
    % ssx = 2^ceil(log2(params(2)));
    % ssy = 2^ceil(log2(params(1)));
    % ss = max(ssx, ssy);
    % W = Wavelet('Daubechies',4,2);
    % % W = Wavelet('Haar',2,3);
    %
    %
    % for n=1:nIter
    %     X = zpad(ifft2c(x),ss,ss,params(3));
    %     X = SoftThresh(W*(X),wavWeight);
    %     x = fft2c(crop(W'*(X),params(2),params(1),params(3))).*(1 - masku) + data;
    % end
    
    toc
    
    save([application(5:end),'/accelfactor_',num2str(accelfactor),'/irls_cs_recon1.mat'], 'kxkykzc','-v7.3');
    
    
    % aba = combinechannels(fft2c(squeeze(kxkykzc(:,:,:,:,15))),4,2); % reconstructed image
    % mip = max(aba(:,:,:),[],3);
    % figure;imshow(abs(mip),[]);
    
end