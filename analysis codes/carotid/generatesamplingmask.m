%% Generates a sampling mask
% cleaning the workspace
clear all; close all; clc;

N = [128, 256];     % imagesize
accelfactor = 8;  % acceleration factor

% pseudo random sampling
pdf = genPDF(N, 5, 1/accelfactor, 2, 0, 0);
mask = genSampling(pdf,500,2);
imshow(mask)
save(['samplingmasks/sampling_mask_',num2str(accelfactor),'.mat'],'mask','pdf');