
clear all; clc, close all

load data_3DsprDE; % contain raw data and bunch of imaging parameters


%% Lienar off-resonance correction turned OFF


f0 = zeros(1,Nz);
for zz=1: Nz
    kcorr(:,:,zz) = k;
    wcorr(:,:,zz) = w;
end

%% SPIRiT recon


% Recon parameters
kSize = [7,7,3];
CalibSize = [ 30,30,8];

nIterCG = 12;                  
CalibTyk = 0.02;               

lambda = 3;                     % Ratio between data and calibration consistency. 1 recommended when density compensation is used.

for zz=1:Nz
        
    FT{zz} = NUFFT(kcorr(:,:,zz),w,[0,0],[ N,N] );

end

sampidx_tot1 = ones((size(kcorr,2)/accel),Nz);
for i = 2:(size(kcorr,2)/accel)
    sampidx_tot1(i,:) = sampidx_tot1(i - 1,:) + accel;
end

for i = 2:accel
    sampidx_tot1(:,i:accel:end) = sampidx_tot1(:,i:accel:end) + (i - 1);
end


FT_SPR3d = SPR3dFULL( FT, [N,N,Nz], [Nk,Nspr,Nz], Nc, sampidx_tot1, f0 );

imund = FT_SPR3d'*(raw.*repmat(sqrt(w),[1,1,Nz,Nc]) );
imundss = sumsquare( permute( imund, [ 4,1,2,3] ));

% Calibration
    
raw_calib = raw;


im_calib = FT_SPR3d'*(raw_calib.*repmat(sqrt(w),[1,1,Nz,Nc]));

clear raw_calib;
kData = fft3c_spirit(im_calib);
clear im_calib;
    
kCalib = crop3d(kData,[CalibSize,Nc]); % crop center k-space for calibration
clear kData;
    
[AtA,] = corrMatrix3d(kCalib,kSize);
   
kernel = zeros([kSize,Nc, Nc]);
for cc=1:Nc
    disp(sprintf('Calibrating coil %d',cc));
    [ kernel_cc rawkernel_cc] = calibrate3d( AtA, kSize, Nc, cc, CalibTyk);
    kernel(:,:,:,:,cc) = kernel_cc;
end
disp('Done Calibrating');
    
% Init SPIRiT operator
GOP = SPIRiT3d(kernel, 'image', [N N Nz] );
    

[ imrec FLAG,RELRES,ITER,RESVEC,LSVEC]= cgNUSPIRiT3d(double(raw).*repmat(sqrt(w),[1,1,Nz,Nc]),imund,FT_SPR3d,GOP,nIterCG,lambda);

imrec_ss = sumsquare( permute( imrec, [ 4,1,2,3] ));

