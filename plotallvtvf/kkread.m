function [k,wt] = kkread(filename,num,numr);

%
% KKREAD	function [k wt] = kkread(filename,num,numr)
%
%	reads a Kevin-King style K-space file.
%	filename is the name of the kspacefile, num is the number
%	of phase encodes, numr is the length of the readout.
%
%	returns two matching matrices, K (complex values), and Weight
%
% ksn 7/2/98
%

fil = fopen(filename,'r','ieee-be');

kx = zeros(numr,num);
ky = zeros(numr,num);
wt = zeros(numr,num);

[dump ct] = fread(fil, [3*numr num], 'float32');

kx = dump(1:3:3*numr,:);
ky = dump(2:3:3*numr,:);
wt = dump(3:3:3*numr,:);

k = kx+(j*ky);

fclose(fil);

%disp(['count = ',num2str(ct)]);

