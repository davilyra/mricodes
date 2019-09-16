function rawdata = rawdataRead(coildata, nSe, nPe, nPeFull, nRo, nAcq, nCoil,extraAcq)
%==============================================================================
% rawdata = rawdata_extractor_2DFBI(nPe, nPeFull, nRO, nAcq, nCoil)
% 
% extract k-space raw data from study file outputs (write SVN_GEN_DATA coildata)
% By Cheng Ouyang, March, 2012, version 1.0

% inputs: (raw data specifications)
% coildata    name of coildata, e.g., proc3_coil0, then coildata = proc3_coil
% nSe:        number of slices
% nPe:        number of acquired pe lines            
% nPeFull     full pe lines, usually 64-128-256-512-...
% nRo         full ro lines, usually 64-128-256-512-...
% nAcq        number of repeats
% nCoil       number of channels
% extraAcq    1 = if there are extra lines acquired, 0 = otherwise

% output: 
% rawdata: k-space data complex, dim = nSe x nPe x nRo x nCoil x nAcq x 2
% 
% usage:  rawdata = rawdata_extractor_2DFBI('proc3_coil',1,88,128,512,17,32);
% proc3_coilx data is in MATLAB/data/2012-03-07-CS-FBI/coildata
%==============================================================================

%%

% raw data specifications
nPeZero = nPeFull - nPe;

%fid = fopen('ECG_Prep_coil0','r');
%aaa = fread(fid,'int32');
%ECG_Prep_coil0 data size (1531904) = nPe (88) x nRo (512) x nAcq (17) x 2(real vs. imaginary)
%data size (1531904) = dataDim x nAcq
%i.e., in study file outputs of coildata, only non-zero/acquired data points (nPe not nPeFull) are stored.

% dataDim =  nSe * nPe * nRo * 2;
realS1 = zeros(nSe, nPeFull, nRo, nCoil, nAcq);
imagS1 = realS1;

fclose('all');
filenameIn = strcat(coildata);
fid = fopen(filenameIn, 'r');
rawS1 = fread(fid,'int32');
fclose(fid);

j = 1;
for se=1:nSe
    for acq = 1:nAcq
        for coil = 1:nCoil            
%             op = sprintf('reading data slice %s coil %s', num2str(se),num2str(coil));
%             disp(op);
            for pe = 1:nPeFull,
                if pe > nPeZero,
                    for ro = 1:nRo,       % read data straight
                        realS1( se, pe, ro, coil, acq ) = rawS1( j   );
                        imagS1( se, pe, ro, coil, acq ) = rawS1( j+1 );
                        j = j+2;
                    end
                end
            end
        end
    end
end

rawdata = complex(realS1, imagS1);
rawdata =squeeze(rawdata);

% figure, imagesc(log10(abs(realS1(:,:,2,1))),[0 log10(5000000)]),
% colormap(gray)
