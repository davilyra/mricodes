clear,clc,close all
 
% % The scan parameters are as follows:
% % 
% % scanner: GE Signa 3T EXCITE HD
% % coil: 4-channel neck array (two elements on each side of the neck)
% % flows of interest: left and right carotid arteries and bifurcations, left and right jugular veins, vertebral arteries, etc)
% % 
% % healthy volunteer
% % heart rate: 105 bpm 
% % spatial resolution: 1.4 x 1.4 x 5 mm
% % number of slices: 5 (contiguous, acquired independently)
% % FOV: 16x16 cm 
% % temporal resolution: 11.7 ms (1 TR)
% % velocity FOV: 240 cm/s (+/- 120 cm/s)
% % velocity resolution: 7.5 cm/s 
% % 32 velocity encodes
% % 43 cardiac phases (90% of cardiac cycle)
% % 8 spiral interleaves
% % variable-density spirals (16~4 cm FOV)

rawpath = './rawdata/'; %directory where the raw data files are located
nslices = 5; %number of slices
ncoils = 4; % number of coils

%%%% load aquisition parameters
filename = sprintf('%sslice1.7',rawpath);
[rawdata,usercv,hdr] = rawloadHD_jfn(filename,[],[],[], 1, [], [], []);
nphases = usercv(4); %number of cardiac phases (i.e., temporal frames)
nVE = usercv(7); %number of velocity encoding steps
nread = usercv(8); %number of readout samples
nintl = usercv(14); %number of spiral interleaves
spiralid = usercv(10); %this number identifies which spiral trajectory was used

%%% load spatial parameters, k-space trajectory, and density compensation weights
[kfile,spatfov,spatres] = getspiralparams(spiralid);
[kxkytraj kxkyweights] = kkread(kfile,nintl,nread);
imagesize = ceil(spatfov/spatres); %image size (in pixels)
kxkyweights=kxkyweights/max(kxkyweights(:)); %normalize the density values
kxkytraj = kxkytraj / (2*max(abs(kxkytraj(:)))); %normalize to [-0.5 , 0.5]
%figure,plot(kxkytraj);axis equal; %plot kspace trajectory

% specify velocity encode, cardiac phase and slice to be reconstructed
velocityencode = 17; %(all 32 velocity encodes are loaded by rawloadHD_jfn, but only one is reconstructed
cardiacphase = 1;
slice = 4;

%specifiy coil elements and spiral interleaves to be reconsctructed
coils = 1:ncoils;
spiralinterleaves = 1:nintl;

%filename for specified slice
filename = sprintf('%sslice%d.7',rawpath,slice); 


%%%%%% USAGE: rawdata = rawloadHD_jfn(filename,[],[],[], 1,coils,phases,spiralinterleaves);
%% coils is a vector with any subset of numbers from 1 to ncoils
%% phases is a vector with any subset of numbers from 1 to nphases
%% spiralinterleaves is a vector with any subset of numbers from 1 to nintl
%%%% rawdata is a matrix with dimension: nread x nVE x length(coils) x length(phases) x length(spiralinterleaves) (in this order)

%load rawdata (this loads data from all velocity encodes)
kxkykvct = permute(... %permute will re-order the matrix dimensions
    rawloadHD_jfn(filename,[],[],[], 1,coils,cardiacphase,spiralinterleaves),... %load data from select coils, cardiac phases and spiral interleaves
    [1 5 2 3 4]); %new matrix dimension order will be: readout, spiral interleave, velocity encode, coil, cardiac phase (hence, kx-ky-kv-c-t)


rawdata = kxkykvct(:,:,nVE/2+1,:,1);

rawdata = permute(rawdata,[1 2 4 3]);

size(rawdata)

ktraj = kxkytraj;
kweights = kxkyweights;

save carotid_sp_8intl_4ch.mat rawdata ktraj kweights imagesize -v6


% %% reconstruct image from a select velocity encode and a select cardiac phase
% xyc = zeros(imagesize,imagesize,ncoils);
% for coil = 4, %loop through the coil elements
%     %load k-space data from a single coil
%     kxkydata = kxkykvct(:,:,velocityencode,coil,cardiacphase);
%     
%     % inverse non-uniform Fourier transform using DrFT (SLOW!!)
%     disp(sprintf('reconstructing image from coil #%d...',coil)),tic
%     xyc(:,:,coil) = idrft(kxkydata(:),kxkyweights(:),kxkytraj(:),imagesize); 
%     disp(sprintf('Done! Elapsed time: %f seconds.',toc))
% end;
% im = combine4channels(xyc(:,:,1),xyc(:,:,2),xyc(:,:,3),xyc(:,:,4)); %combine images from 4 channels
% 
% %display result
% normim = flipud(abs(im).'); % should we apply fliplr too?
% normim = normim/max(normim(:));
% imshow(normim)