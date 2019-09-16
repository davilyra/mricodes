clear,clc,close all
 
rawpath = './rawdata/'; %directory where the raw data files are located
nslices = 5; %number of slices
ncoils = 4; % number of coils

%%%% load aquisition parameters
filename = sprintf('%sslice1.7',rawpath);
[rawdata,usercv,hdr] = rawloadHD_jfn(filename,[],[],[], 1, [], [], []);
nphases = usercv(4); %number of cardiac phases (i.e., temporal frames)
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

for slice=1:nslices,
    
    disp(sprintf('converting data from slice #%d...',slice))

	%filename for specified slice
	filename = sprintf('%sslice%d.7',rawpath,slice); 
	
	%%%%%% USAGE: rawdata = rawloadHD_jfn(filename,[],[],[], 1,coils,phases,spiralinterleaves);
	%% coils is a vector with any subset of numbers from 1 to ncoils
	%% phases is a vector with any subset of numbers from 1 to nphases
	%% spiralinterleaves is a vector with any subset of numbers from 1 to nintl
	%%%% rawdata is a matrix with dimension: nread x nVE x length(coils) x length(phases) x length(spiralinterleaves) (in this order)
	
	%loads and reorders entire dataset
	kxkykvct = permute(... %permute will re-order the matrix dimensions
        rawloadHD_jfn(filename,[],[],[], 1,1:ncoils,1:nphases,1:nintl),...
        [1 5 2 3 4]); %new matrix dimension order will be: readout, spiral interleave, velocity encode, coil, cardiac phase (hence, kx-ky-kv-c-t)(
    
    filename = sprintf('%sslice%d.mat',rawpath,slice);
    save(filename,'kxkykvct','kxkytraj','kxkyweights','imagesize')

end;

disp('All done!')