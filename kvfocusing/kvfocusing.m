function [uddata,udkv] = kvfocusing(vddata,vdkv,varwidthsinc,dynFOVcenter,usehamm,overgridfactor)
% Kv-FOCUSING: FOV Centering Using Sinc INterpolation along Kv
%
% Joao L. A. Carvalho, 2007 - ISMRM 2007 abstract #2514
% 
% INPUTS:
% vddata: FVE data
% vdkv: Kv locations of the FVE data samples (e.g., in s/cm)
% varwidthsinc: set to 1 to use variable-width sinc kernels, 0 for constant-width. Recommended: 1.
% dynFOVcenter: set to 1 to use dynamic FOV centering, 0 otherwise. Recommended: 1.
% usehamm: set to 1 to use Hamming windowing, 0 otherwise. Recommended: 1.
% overgridfactor: over-grid factor. Recommended: 1.5.
%
% OUTPUTS:
% uddata: interpolated uniform-density FVE data
% udkv: uniform-density Kv locations

% sorts the data in ascending order of kv-location
[vdkv,newkvorder] = sort(vdkv);
vddata = vddata(newkvorder);

%calculates minimum sample-interval
mindeltakv = min(diff(vdkv));
if mindeltakv == 0,
    error('only 1 sample per kv location');
end;

%computes weighting function
weight = voronoi1d(vdkv); 

if(dynFOVcenter),
    %centers the FOV based on a phase-contrast reference
    pcref = calcpcref(vdkv,vddata);
else
    %no FOV recentering
    pcref = 0;
end;

%builds uniform sampling grid
maxkv = max(abs(vdkv));
mindeltakv = maxkv/ceil(maxkv/mindeltakv);
% temp = -(mindeltakv:mindeltakv:(maxkv+mindeltakv))'; %negative side
% temp = temp(end:-1:1); %reverse order
% udkv = [temp;(0:mindeltakv:maxkv)']; %concatenates with positive side
udkv = (-maxkv):mindeltakv:(maxkv-mindeltakv);
nVE = length(udkv);
maxkv = -udkv(1);            
ktraj = vdkv/maxkv/2;

%resamples data onto new grid
[uddata,udkv] = gridsinc1d(ktraj,vddata,weight,overgridfactor,nVE,varwidthsinc,pcref,maxkv,usehamm);
udkv = udkv/nVE*maxkv*2;
