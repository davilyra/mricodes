addpath(genpath('/home/mri-vh/Davi/functions'))
addpath(genpath('/home/mri-vh/Davi/ADMM_CS_SENSE'))
addpath(genpath('/home/mri-vh/Davi/Marc_SPIRiT'))

clear all;
close all;
clc;

load('params.mat');

type = 1; % [1] = pese; [2] = pero

accelfactors = [2];
% accelfactors = [2 4 5 8 10];
uslth = length(accelfactors);

for iii = 1:uslth
    accelfactor = accelfactors(iii);
    
    switch type
        case 1
            application = 'FBI/pese';
            mkdir([application(5:end),'/accelfactor_',num2str(accelfactor)]);
            load([application(5:end),'/accelfactor_',num2str(accelfactor),...
                '/rawdata_fbi1_',application(5:end),'_us_',num2str(accelfactor),'.mat']);
            
            load(['/home/mri-vh/Davi/samplingmasks/',application,'/sampling_mask_',num2str(accelfactor),'.mat']);
            mask = permute(repmat(mask,[1,1,params(2),params(4),params(5)]),[3 1 2 4 5]);
            
            N = params(1);
            
            N2 = params(1)*params(3);
            NL = params(1)*params(3)*params(4);
            
            % Total variance in two direction, a function version is also available
            T1 = sparse([1:NL,2:NL],[1:NL,1:NL-1],[ones(1,NL),-ones(1,NL-1)]); % construct T1
            
            T2 = sparse([1:NL,N2 + 1:NL],[1:NL,1:N2*(params(4)-1)],...
                [-ones(1,NL),ones(1,N2*(params(4)-1))]); % construct T2
            
            D = [T1;T2]; % construct tall matrix D
            
            gama = 0.001; % dummy penalty, needs improvement
            beta = 0.005;
            opt.kx = params(1);
            opt.ky = params(3);
            opt.coils = 16;
            opt.itermax = 10;
            opt.tol = 1e-4;
            opt.gama = 0.1;
            
%             wavWeight = 0.03125;
%             nIter = 15;
            
            kxkykzc = zeros(params(2),params(1),params(3),params(4),params(5));
            for bbti = 15;
%             for bbti = 1:params(5);
                for ii = 1:params(2)
                    kU = squeeze(rd_fbi_us(ii,:,:,:,bbti));
                    masku = squeeze(mask(ii,:,:,:,bbti));
                    
                    b = gama*(ifft2c(masku.*kU)); % gama*DFT'*U'*kU to prepare for A'*b
                    b = b(:);
                    
                    % set initial values
                    xk = b;
                    xk_p = xk;
                    xk1 = zeros(NL,1);
                    vk = zeros(2*NL,1); % dummy variable vk,uk
                    iter = 1;
                    itermax = 20;
                    tol = 1e-4;
                    
                    tic,
                    while iter<itermax && norm(xk1-xk_p)>tol
                        dk1 = shrink1(vk + D*xk, beta*gama);
                        uk1 = vk + D*xk - dk1;
                        [xk1,count] = CS_SENSE_recon(ones(1,1,opt.coils), masku, D, b + D'*(dk1 - uk1), xk, opt);
                        
                        xk_p = xk;
                        xk = xk1;
                        dk = dk1;
                        vk = uk1;
                        
                        %                 imagesc(combinechannels(reshape(abs(xk),opt.kx,opt.ky,opt.coils),3,2)); colormap('gray');xlabel(num2str(ii)); drawnow
                        iter = iter + 1;
                    end
                    
                    kxkykzc(ii,:,:,:,bbti) = reshape(xk1,opt.kx,opt.ky,opt.coils);
                    
                    %             ssx = 2^ceil(log2(opt.kx));
                    %             ssy = 2^ceil(log2(opt.ky));
                    %             ss = max(ssx, ssy);
                    %             W = Wavelet('Daubechies',4,2);
                    
                    x = ifft2c(squeeze(kxkykzc(ii,:,:,:,bbti)));
                    
                    %             data = x;
                    %
                    %             for n=1:nIter
                    %                 % apply wavelet thresholding
                    %                 X = zpad(ifft2c(x),ss,ss,opt.coils);
                    %                 X = SoftThresh(W*(X),wavWeight);
                    %                 x = fft2c(crop(W'*(X),opt.kx,opt.ky,opt.coils)).*(1 - masku) + data;
                    %
                    %             end
                    %
                    iter = 1;
                    
                    %             kxkykzc(ii,:,:,:) = fft2c(x);
                    kxkykzc(ii,:,:,:,bbti) = (x);
                    
                    toc,
                end
                
                aba = combinechannels(fft2c(squeeze(kxkykzc(:,:,:,:,bbti))),4,2); % reconstructed image
                mip = max(aba(:,:,:),[],3);
                figure;imshow(abs(mip),[]);
                
            end
%             save([application(5:end),'/accelfactor_',num2str(accelfactor),'/admm_cs_recon1.mat'], 'kxkykzc','-v7.3');            
            
            
        case 2
            application = 'FBI/pero';
            mkdir([application(5:end),'/accelfactor_',num2str(accelfactor)]);
            load('params.mat');
            
            load([application(5:end),'/accelfactor_',num2str(accelfactor),...
                '/rawdata_fbi_',application(5:end),'_us_',num2str(accelfactor),'.mat']);
            
            
            load(['/home/mri-vh/Davi/samplingmasks/',application,'/sampling_mask_',num2str(accelfactor),'.mat']);
            mask = permute(repmat(mask,[1,1,params(3),params(4),params(5)]),[2 1 3 4 5]);
            
            N = params(1);
            
            N2 = params(1)*params(2);
            NL = params(1)*params(2)*params(4);
            
            
            % Total variance in two direction, a function version is also available
            T1 = sparse([1:NL,2:NL],[1:NL,1:NL-1],[ones(1,NL),-ones(1,NL-1)]); % construct T1
            
            T2 = sparse([1:NL,N2 + 1:NL],[1:NL,1:N2*(params(4)-1)],...
                [-ones(1,NL),ones(1,N2*(params(4)-1))]); % construct T2
            
            D = [T1;T2]; % construct tall matrix D
            
            gama = 0.001; % dummy penalty, needs improvement
            beta = 0.005;
            opt.kx = params(2);
            opt.ky = params(1);
            opt.coils = 16;
            opt.itermax = 10;
            opt.tol = 1e-4;
            opt.gama = 0.1;
            
            %         wavWeight = 0.03125;
            %         nIter = 15;
            
            kxkykzc = zeros(params(2),params(1),params(3),params(4),params(5));
            for bbti = 1:params(5);
                for ii = 1:params(3)
                    kU = squeeze(rd_fbi_us(:,:,ii,:,bbti));
                    masku = squeeze(mask(:,:,ii,:,bbti));
                    
                    %             img1 = squeeze(combinechannels(fft2c(squeeze(rd_fbi_us(:,:,ii,:,bbti))),3,2));
                    
                    b = gama*(ifft2c(masku.*kU)); % gama*DFT'*U'*kU to prepare for A'*b
                    b = b(:);
                    
                    % set initial values
                    xk = b;
                    xk_p = xk;
                    xk1 = zeros(NL,1);
                    vk = zeros(2*NL,1); % dummy variable vk,uk
                    iter = 1;
                    itermax = 20;
                    tol = 1e-4;
                    
                    tic,
                    while iter<itermax && norm(xk1-xk_p)>tol
                        dk1 = shrink1(vk + D*xk, beta*gama);
                        uk1 = vk + D*xk - dk1;
                        [xk1,count] = CS_SENSE_recon(ones(1,1,opt.coils), masku, D, b + D'*(dk1 - uk1), xk, opt);
                        
                        % this version don't use the FT class
                        
                        xk_p = xk;
                        xk = xk1;
                        dk = dk1;
                        vk = uk1;
                        
                        iter = iter + 1;
                    end
                    
                    kxkykzc(:,:,ii,:,bbti) = reshape(xk1,opt.kx,opt.ky,opt.coils);
                    
                    %             ssx = 2^ceil(log2(opt.kx));
                    %             ssy = 2^ceil(log2(opt.ky));
                    %             ss = max(ssx, ssy);
                    %             W = Wavelet('Daubechies',4,2);
                    
                    x = ifft2c(squeeze(kxkykzc(:,:,ii,:,bbti)));
                    %     imshow(abs(x(:,:,1)),[]), pause(.5)
                    % %             data = x;
                    %
                    %             for n=1:nIter
                    %                 % apply wavelet thresholding
                    %                 X = zpad(ifft2c(x),ss,ss,opt.coils);
                    %                 X = SoftThresh(W*(X),wavWeight);
                    %                 x = fft2c(crop(W'*(X),opt.kx,opt.ky,opt.coils)).*(1 - masku) + data;
                    %
                    %             end
                    
                    iter = 1;
                    
                    kxkykzc(:,:,ii,:,bbti) = (x);
                    
                    toc,
                    
                end
            end
            save([application(5:end),'/accelfactor_',num2str(accelfactor),'/admm_cs_recon.mat'], 'kxkykzc','-v7.3');
            %         aba = combinechannels(fft2c(squeeze(kxkykzc(:,:,:,:))),4,2); % reconstructed image
            %         mip = max(aba(:,:,:),[],3);
            %         figure;imshow(abs(mip),[]);
            
    end
    
end