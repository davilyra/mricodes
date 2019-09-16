function displayvtvf(vt,optr,maxveloc,dbrange,X,Y,slice,displayvfspace)

nVE = size(vt,1);
nphases = size(vt,2);

taxis = (0:(nphases-1))*optr/1000;
vaxis = (((0:(nVE-1))/nVE)-0.5)*maxveloc*2;
faxis = ((0:(nphases-1))/nphases*(1000000/optr))-(1000000/optr/2);

%plots the time-velocity distribution
subplot(223),
imagesc(taxis,vaxis,abs(vt))
colormap(gray)
set(gca,'YDir','normal')
xlabel('time (ms)')
ylabel('velocity (cm/s)')
title(sprintf('time velocity distribution for (x=%d,y=%d,z=%d)',round(X),round(Y),slice))

if displayvfspace,
    %plots v-f space
    vf = fftshift(abs(fft(vt,[],2)),2);
    subplot(224),
    if dbrange==0,
        imagesc(faxis,vaxis,vf/max(vf(:)),[0 1]),
    else
        imagesc(faxis,vaxis,20*log10(vf/max(vf(:))),[-dbrange 0]),
    end;
    set(gca,'YDir','normal'),
    colorbar
    xlabel('frequency (Hz)')
    ylabel('velocity (cm/s)')
end;