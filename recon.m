clear,clc,close all

load carotid_sp_8intl_4ch.mat %rawdata ktraj kweights imagesize

[nread,nintl,ncoils] = size(rawdata);

% rawdata = rawdata(:,1:2:nintl,:);
% ktraj = ktraj(:,1:2:nintl);
% kweights = kweights(:,1:2:nintl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STARTING UP THE NUFFT ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd fessler; setup; cd ..
Nd = [imagesize,imagesize];
Jd = [6,6];
overgridfactor = 2;
om(:,1) = 2*pi*real(ktraj(:)); %kspace trajectory
om(:,2) = 2*pi*imag(ktraj(:));
st = nufft_init(om, Nd, Jd, overgridfactor*Nd, Nd/2); %initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xyc = zeros(imagesize,imagesize,ncoils);
for coil = 1:ncoils,    
    weighteddata = rawdata(:,:,coil).*kweights; %weights the kspace data based on sampling density
    m = nufft_adj(weighteddata(:),st)/imagesize; %non-cartesian inverse Fourier transform along kx-ky      
    xyc(:,:,coil) = m.'; %flips the image and stores it in xyc
end;   

%combines the images from the 4 coils into a single image
m = combine4channels(xyc(:,:,1),xyc(:,:,2),xyc(:,:,3),xyc(:,:,4));
m = abs(m);
m = m/max(m(:));

imshow(m)
set(gca,'YDir','normal')
