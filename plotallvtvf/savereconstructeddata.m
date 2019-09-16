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


%% loading data
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

% analysis points
axesX = [74.8667 77.9333 35 33.4667 28.8667];
axesY = [71.0333 66.9444 75.1222 70.5222 65.9222];

% saving data (SPIRiT)
xyspisave = abs(xyspi)/(max(abs(xyspi(:))));
imwrite(abs(xyspisave),['xyspi_',t,'_accel_',ta,'.png'],'PNG');

% saving data (nuFFT SOS)
xysave = abs(xy)/(max(abs(xy(:))));
imwrite(abs(xysave),['xynufft_',t,'_accel_',ta,'.png'],'PNG');


%% takes the velocity distribution from the prescribed pixel
% SPIRiT
vtspi_1 = xyvtspi(round(axesX(1)),round(axesY(1)),:,:); % Left Carotid Bifurcation
vtspi_1 = permute(vtspi_1,[3 4 1 2]);

vtspi_2 = xyvtspi(round(axesX(2)),round(axesY(2)),:,:); % Left Jugular Vein
vtspi_2 = permute(vtspi_2,[3 4 1 2]);

vtspi_3 = xyvtspi(round(axesX(3)),round(axesY(3)),:,:); % Right External Carotid Artery
vtspi_3 = permute(vtspi_3,[3 4 1 2]);

vtspi_4 = xyvtspi(round(axesX(4)),round(axesY(4)),:,:); % Right Internal Carotid Artery
vtspi_4 = permute(vtspi_4,[3 4 1 2]);

vtspi_5 = xyvtspi(round(axesX(5)),round(axesY(5)),:,:); % Right Jugular Vein
vtspi_5 = permute(vtspi_5,[3 4 1 2]);

% nuFFT SOS
vt_1 = xyvt(round(axesX(1)),round(axesY(1)),:,:); % Left Carotid Bifurcation
vt_1 = permute(vt_1,[3 4 1 2]);

vt_2 = xyvt(round(axesX(2)),round(axesY(2)),:,:); % Left Jugular Vein
vt_2 = permute(vt_2,[3 4 1 2]);

vt_3 = xyvt(round(axesX(3)),round(axesY(3)),:,:); % Right External Carotid Artery
vt_3 = permute(vt_3,[3 4 1 2]);

vt_4 = xyvt(round(axesX(4)),round(axesY(4)),:,:); % Right Internal Carotid Artery
vt_4 = permute(vt_4,[3 4 1 2]);

vt_5 = xyvt(round(axesX(5)),round(axesY(5)),:,:); % Right Jugular Vein
vt_5 = permute(vt_5,[3 4 1 2]);


%% errors
% SPIRiT and nuFFT
xyspien = abs(xyspi)/energyd(abs(xyspi));
xyen = abs(xy)/energyd(abs(xy));

error0 = abs(xyen - xyspien);

% flow profiles
vt_1en = abs(vt_1)/energyd(abs(vt_1));
vt_2en = abs(vt_2)/energyd(abs(vt_2));
vt_3en = abs(vt_3)/energyd(abs(vt_3));
vt_4en = abs(vt_4)/energyd(abs(vt_4));
vt_5en = abs(vt_5)/energyd(abs(vt_5));

vtspi_1en = abs(vtspi_1)/energyd(abs(vtspi_1));
vtspi_2en = abs(vtspi_2)/energyd(abs(vtspi_2));
vtspi_3en = abs(vtspi_3)/energyd(abs(vtspi_3));
vtspi_4en = abs(vtspi_4)/energyd(abs(vtspi_4));
vtspi_5en = abs(vtspi_5)/energyd(abs(vtspi_5));

error1 = abs(vt_1en - vtspi_1en);
error2 = abs(vt_2en - vtspi_2en);
error3 = abs(vt_3en - vtspi_3en);
error4 = abs(vt_4en - vtspi_4en);
error5 = abs(vt_5en - vtspi_5en);


%%  saving into images
imwrite((flipud(abs(vtspi_1)/max(abs(vtspi_1(:))))),['vtspi_',t,'_',ta,'_LCB.png'],'PNG');
imwrite((flipud(abs(vtspi_2)/max(abs(vtspi_2(:))))),['vtspi_',t,'_',ta,'_LJV.png'],'PNG');
imwrite((flipud(abs(vtspi_3)/max(abs(vtspi_3(:))))),['vtspi_',t,'_',ta,'_RECA.png'],'PNG');
imwrite((flipud(abs(vtspi_4)/max(abs(vtspi_4(:))))),['vtspi_',t,'_',ta,'_RICA.png'],'PNG');
imwrite((flipud(abs(vtspi_5)/max(abs(vtspi_5(:))))),['vtspi_',t,'_',ta,'_RJV.png'],'PNG');

%  saving into image
imwrite((flipud(abs(vt_1)/max(abs(vt_1(:))))),['vtnufft_',t,'_',ta,'_LCB.png'],'PNG');
imwrite((flipud(abs(vt_2)/max(abs(vt_2(:))))),['vtnufft_',t,'_',ta,'_LJV.png'],'PNG');
imwrite((flipud(abs(vt_3)/max(abs(vt_3(:))))),['vtnufft_',t,'_',ta,'_RECA.png'],'PNG');
imwrite((flipud(abs(vt_4)/max(abs(vt_4(:))))),['vtnufft_',t,'_',ta,'_RICA.png'],'PNG');
imwrite((flipud(abs(vt_5)/max(abs(vt_5(:))))),['vtnufft_',t,'_',ta,'_RJV.png'],'PNG');

imwrite((flipud(abs(error0)/max(abs(error0(:))))),['error_spirit_nufft_',t,'_us_',ta,'.png'],'PNG');
imwrite((flipud(abs(error1)/max(abs(error1(:))))),['error_spirit_nufft_',t,'_us_',ta,'_LCB.png'],'PNG');
imwrite((flipud(abs(error2)/max(abs(error2(:))))),['error_spirit_nufft_',t,'_us_',ta,'_LJV.png'],'PNG');
imwrite((flipud(abs(error3)/max(abs(error3(:))))),['error_spirit_nufft_',t,'_us_',ta,'_RECA.png'],'PNG');
imwrite((flipud(abs(error4)/max(abs(error4(:))))),['error_spirit_nufft_',t,'_us_',ta,'_RICA.png'],'PNG');
imwrite((flipud(abs(error5)/max(abs(error5(:))))),['error_spirit_nufft_',t,'_us_',ta,'_RJV.png'],'PNG');

% max(error0(:)), max(error1(:)), max(error2(:)), max(error3(:)), max(error4(:)), max(error5(:))