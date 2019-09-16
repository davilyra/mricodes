function  res = SPR3dFULL(FT_NC, imSize, dataSize, nCoil, sampidx, f0 )
% res = SPR3dFULL(FT_NC, imSize, dataSize, nCoil, sampidx, f0 )
% Inputs:
%   - FT_NC: NUFFT Object
%   - imSize: Output image size
%   - dataSize: Input k-space data size, e.g.: [# phase encoding lines, #
%   freq encoding lines, # levels in kz], [# spiral readouts, # spiral
%   interleaves, # velocity encoding planes].
%   - nCoil: number of coils used in the data acquisition
%   - sampidx: a kind of sampling matrix, for spiral readouts it has the
%   followin size: (# of interleaves effectively used,
%   # of levels in kz or kv). Example:
%       sampidx = ones((nintl/undersamplingfactor),nVE);
%       for i = 2:(nintl/undersamplingfactor)
%           sampidx(i,:) = sampidx(i - 1,:) + undersamplingfactor;
%       end
%       for i = 2:undersamplingfactor
%           sampidx(:,i:undersamplingfactor:end) = sampidx(:,i:undersamplingfactor:end) + (i - 1);
%       end
%   (obtained from Taehoon Shin's original example)
%
%   - f0: vector with zeros, size = [1, # kv's or kz's]
%
% Outputs:
%   - res: SPIRiT 3D Object

% Taehoon Shin 2011
% Edited by Davi M. Lyra-Leite, Apr 2013.

res.adjoint = 0;
% res.kmask = mask;
% res.imSize = imSize;
% res.dataSize = size(mask);
res.FT = FT_NC;
res.imSize = imSize;
res.dataSize = dataSize;
res.Nc = nCoil;
res.sampidx = sampidx;
res.f0 = f0;
res = class(res,'SPR3dFULL');
