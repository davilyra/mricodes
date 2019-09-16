clear,clc,close all

rawpath = './rawdata/';


nphases = 43;
ncoils = 4;
nintl = 8;

coils = 1:ncoils;
phases = 1:nphases;

%for slice=1:5,
slice = 1;

filename = sprintf('%sslice%d.7',rawpath,slice);
raw = rawloadHD_jfn(filename,[],[],[], 1,coils,phases,1:nintl);
%  1012          32           4          43           8
% readout        Kv         coils       time        interleaves

