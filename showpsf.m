function showpsf(xyvt,maxveloc,optr,spatfov,dbrange)

nVE = size(xyvt,3);
nphases = size(xyvt,4);

xyvf = abs(fftshift(fft(xyvt,[],4),4));
maxxyvf = max(xyvf(:));

figure(1),
for v=1:nVE,
    for p=1:nphases,
        subplot(nVE,nphases,p+(v-1)*nphases),
        if dbrange==0,
            imagesc([],[],xyvf(:,:,v,p)/maxxyvf,[0 1])
        else
            imagesc([],[],20*log10(xyvf(:,:,v,p)/maxxyvf),[-dbrange 0]),
        end;
        set(gca,'YDir','normal'),
        axis square
        %colorbar
        f = ((p-1)/nphases-0.5)*(1000000/optr);
        vv = ((v-1)/nVE-0.5)*maxveloc*2;
        set(gca,'XTick',[],'YTick',[])
        if(v==1),title(sprintf('f = %0.2f Hz',f));end;
        if(p==1),ylabel(sprintf('v = %0.0f cm/s',vv)); end;
        colormap(gray)
    end;
end;

figure(2),
nx = size(xyvf,1);
dx = spatfov/nx;
x = (-spatfov/2):dx:(spatfov/2-dx);
xcentral = floor(nx/2+1);
for v=1:nVE,
    for p=1:nphases,
        subplot(nVE,nphases,p+(v-1)*nphases),
        if dbrange==0,
            plot(x,xyvf(xcentral,:,v,p)/maxxyvf),
            axis([-spatfov/2 spatfov/2 0 1])
            %set(gca,'YTick',0:.2:1)
        else
            plot(x,20*log10(xyvf(xcentral,:,v,p)/maxxyvf)),
            axis([-spatfov/2 spatfov/2 -dbrange 0])
            %set(gca,'YTick',-dbrange:10:0)
        end;
        f = ((p-1)/nphases-0.5)*(1000000/optr);
        vv = ((v-1)/nVE-0.5)*maxveloc*2;
        if(v==1),title(sprintf('f = %0.2f Hz',f));end;
        if(p==1),ylabel(sprintf('v = %0.0f cm/s',vv)); end;
        grid
    end;
end;