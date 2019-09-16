function [filename,maxveloc,optr,nphases,nVE,nread,nintl,spiralid,spatfov,spatres,pixels,kxkytraj,kxkyweights]=readdataparams(rawpath,slice)

%this loads the usercv variables from the rawdata file
filename = sprintf('%sslice%d.7',rawpath,slice); %filename of the current slice
[rawdata,usercv,hdr] = rawloadHD_jfn(filename,[],[],[], 1, [], [], []);
maxveloc = usercv(1); % maximum velocity value (1/2 of velocity field-of-view)
optr = usercv(2); % TR duration (in microseconds)
nphases = usercv(4); % number of cardiac phases (i.e., temporal frames)
nVE = usercv(7); % number of velocity encoding steps
nread = usercv(8); % number of readout samples
nintl = usercv(14); % number of spiral interleaves
spiralid = usercv(10); % this number identifies which spiral trajectory was used
% heartrate = usercv(5); % heart rate of the subject during the scan
% RRpct = usercv(6); % percentage of cardiac cycle that was covered
% flipangle = usercv(15); % flip angle (in degrees)
% vesperbeat = usercv(3);
% vres = usercv(9);
% vegrad = usercv(11);
% realtime = usercv(12);
% nocine = usercv(13);
% ia_rf1 = usercv(16); 
% densityreductionfactor = usercv(17);
% oblique = usercv(18);
% intlsperbeat = usercv(19);

%==========================================================================
% load spatial parameters
switch(spiralid)
    case 19,
        kfile='recon16cm14mm8intl4vd';
        spatfov = 160;
        spatres = 1.4;         
    otherwise,
        error(['unexpected spiralid value: ',num2str(spiralid)]);
end;
[kxkytraj kxkyweights] = kkread(kfile,nintl,nread);

%nread = 1012;
kxkytraj = kxkytraj(1:nread,:);
kxkyweights = kxkyweights(1:nread,:);
kxkytraj = kxkytraj*256/spatfov; %corrects k-space trajectory
spatres = 1/(2*abs(kxkytraj(end,1))); %calculates spatial resolution from kspace trajectory
%figure,plot(kxkytraj);axis equal; %plots kspace trajectory
kxkytraj = kxkytraj*spatres; %normalize to [-0.5 , 0.5]
kxkyweights = kxkyweights/max(kxkyweights(:));%/sqrt(2); %normalizes the density values to [0 , sqrt(2)/2]

pixels = ceil(spatfov/spatres);      % image size