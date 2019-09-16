function initialization

%%spirit and nufft initialization
if exist('nufftinit.m')==0
	addpath('.\SPIRiT_v0.1\utils')
    addpath('SPIRiT_v0.1')
    cd fessler,setup, cd ..
    %addpath('.\SPIRiT_v0.1\nufft_files');
end
% if exist('nufft.m')==0
%     addpath('.\SPIRiT_v0.1\nufft_files');
% end