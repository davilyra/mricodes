% RECONSTRUCTION OF THE TIME-VELOCITY DISTRIBUTION FROM
% MULTIPLE VOXELS IN A SLICE, WHICH ARE MANUALLY
% PRESCRIBED BY THE USER. CLOSE THE FIGURE OR PRESS CTRL-C TO STOP.
%
%
% Written by Joao L. A. Carvalho <joaoluiz@gmail.com>
% Department of Electrical Engineering
% University of Brasilia, Brazil
%
% July 14, 2008
%
% Changes by Davi Marco Lyra Leite <davi@ieee.org>
% November 3rd, 2011


clear,clc,close all

disp('loading data...')
t = input('Com qual slice você deseja trabalhar (1 a 5)? \n R: ','s');
t1 = str2double(t);

ta = input('Com qual taxa de subamostragem você deseja trabalhar (1, 2 ou 4)? \n R: ','s');
t2 = str2double(ta);

% for SPIRiT data
cd datapaper\
load(['slice_',t,'_',ta,'.mat'])
load(['slice_',t,'_',ta,'_xy.mat'])

xyvtspi = xyvt;
maxvelocspi = maxveloc;
optrspi = optr;
xyspi = xy;

% for nuFFT SOS data
load(['slice_',t,'_sosnufft_',ta,'.mat'])
cd ..


% plotting data (SPIRiT)
nVE = size(xyvtspi,4);
nphases = size(xyvtspi,5);

disp('Plotting data using SPIRiT. Close the figure or press Ctrl-C to stop.')
figure,

while(1),
    
    %plots the magnitude image (central kv, first cardiac phase)
    subplot(211),
    imshow(abs(xyspi),[ ])
    set(gca,'YDir','normal')
    title(sprintf('slice %d and undersamplingfactor %d',t1,t2))

    %user will click on a pixel of the image
    xlabel('CLICK ON A BLOOD VESSEL / CLOSE FIGURE TO STOP')
    [Y,X] = ginput(1); %X, Y,

    %takes the velocity distribution from the prescribed pixel
    vtspi = xyvtspi(round(X),round(Y),:,:);
    vtspi = permute(vtspi,[3 4 1 2]);

    %plots the time-velocity distribution
    taxis = (0:(nphases - 1))*optrspi/1000;
    vaxis = (((0:(nVE-1))/nVE)-0.5)*maxvelocspi*2;
    subplot(212),imagesc(taxis,vaxis,abs(vtspi))
    colormap(gray)
    set(gca,'YDir','normal')
    xlabel('time (ms)')
    ylabel('velocity (cm/s)')

end;


% plotting data (nuFFT)
% nVE1 = size(xyvt,4);
% nphases1 = size(xyvt,5);
% 
% disp('Plotting data using nuFFT. Close the figure or press Ctrl-C to stop.')
% figure,
% 
% while(1),
%     
%     %plots the magnitude image (central kv, first cardiac phase)
%     subplot(211),
%     imshow(abs(xy),[ ])
%     set(gca,'YDir','normal')
%     title(sprintf('slice %d and undersamplingfactor %d',t1,t2))
% 
%     %user will click on a pixel of the image
%     xlabel('CLICK ON A BLOOD VESSEL / CLOSE FIGURE TO STOP')
%     [Y,X] = ginput(1); X, Y,
% 
%     %takes the velocity distribution from the prescribed pixel
%     vt = xyvt(round(X),round(Y),:,:);
%     vt = permute(vt,[3 4 1 2]);
% 
%     %plots the time-velocity distribution
%     taxis = (0:(nphases1 - 1))*optr/1000;
%     vaxis = (((0:(nVE1 - 1))/nVE1) - 0.5)*maxveloc*2;
%     subplot(212),imagesc(taxis,vaxis,abs(vt))
%     colormap(gray)
%     set(gca,'YDir','normal')
%     xlabel('time (ms)')
%     ylabel('velocity (cm/s)')
% 
% end;