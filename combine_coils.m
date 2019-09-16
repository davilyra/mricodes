function combineddata = combine_coils(inputdata,coils2recon)
% function combineddata = combine_coils(inputdata,coils2recon)
%
% Inputs:
%   - inputdata: tensor containing as its first dimensions the [x,y] values
%   of the image, third dimension contains either the kv or the v data,
%   fourth dimension is the coils
%   - coils2recon: which coils over the four the user wants to combine
%
% Outputs:
%   - combineddata: tensor of size [x,y,kv or v] that contains the combined
%   information from the coils

ncoils = length(coils2recon);

switch ncoils  % combines the data from the desired coils
    case 1
        combineddata = combine4channels(inputdata(:,:,:,coils2recon),0,0,0);
    case 2
        combineddata = combine4channels(inputdata(:,:,:,coils2recon(1)),inputdata(:,:,:,coils2recon(2)),0,0);
    case 3
        combineddata = combine4channels(inputdata(:,:,:,coils2recon(1)),inputdata(:,:,:,coils2recon(2)),inputdata(:,:,:,coils2recon(3)),0);
    otherwise
        combineddata = combine4channels(inputdata(:,:,:,1),inputdata(:,:,:,2),inputdata(:,:,:,3),inputdata(:,:,:,4));
end