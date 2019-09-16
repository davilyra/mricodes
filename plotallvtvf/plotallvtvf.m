clear,clc,close all

dbrange = 40; %dynamic range (in db) for displaying the v-f space (recommended value: 40). Use dbrange = 0 for linear scale.
slices2recon = [4]; %examples: [4] or [1 2 3 4 5]
voxel2recon = [74,71]; %[X,Y], or [ ] to select from image in while loop (CCA: [74,71])
accelfactors = [4]; %[1 2 4];
vieworderings = [0 1 11 21 2 12 22 3 13 23 4];

X = round(voxel2recon(1));
Y = round(voxel2recon(2));

nslices = length(slices2recon);
nvieworderings = length(vieworderings);
for s = 1:nslices,
    slice = slices2recon(s);
    for vw = 1:nvieworderings,
        viewordering = vieworderings(vw);
        for us = 1:length(accelfactors),
            undersamplingfactor = accelfactors(us);
            load(sprintf('vt_slice%d_x%d_y%d_vw%d_usf%d.mat',slice,X,Y,viewordering,undersamplingfactor),'xy','vt','maxveloc','optr');
            figure,subplot(211),
            imshow(abs(xy),[ ])
            set(gca,'YDir','normal')
            title(sprintf('slice %d, viewordering %d,undersamplingfactor %d',slice,viewordering,undersamplingfactor))
            displayvtvf(vt,optr,maxveloc,dbrange,X,Y,slice,1);
        end;
    end;
end;

