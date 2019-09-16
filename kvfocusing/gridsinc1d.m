function [gridData,k] = gridsinc1d(ktraj,data,weight,overgridfactor,samples,varkernel,pcref,maxkv,usehamm)
% Kv-FOCUSING: FOV Centering Using Sinc INterpolation along Kv
%
% Joao L. A. Carvalho, 2007 - ISMRM 2007 abstract #2514

klocation = samples*ktraj;
data = data.*weight;

gridsize = ceil(samples*overgridfactor);
deltak = samples/gridsize;
k = ((-samples/2):deltak:(samples/2-deltak))';

lenData = length(data);
gridData = zeros(size(k));

convwidth = 1;

udkv = k/samples*maxkv*2;

%disp('Gridding data...');
for p=1:lenData,
    %disp(sprintf('%d/%d',p,lenData));
    if usehamm==1, myhamwin = 0.53836 - 0.46164*cos(2*pi*(k-klocation(p)-samples/2)/samples);
    else myhamwin = 1;
    end;
    if(varkernel),
        convwidth = weight(p); 
    end;
    kernel = sinc((k-klocation(p))/convwidth);
    kvlocp = klocation(p)/samples*maxkv*2;
    kernel = kernel.*exp(-j*2*pi*(udkv-kvlocp)*pcref);
    gridData = gridData + data(p)*1/convwidth*kernel.*myhamwin;
end;