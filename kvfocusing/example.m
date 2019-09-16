load vdfvedata.mat

overgridfactor = 1.5;

%%%% KVFOCUSING
[uddata,udkv] = kvfocusing(vddata,vdkv,1,1,1,overgridfactor);
veldistr = abs(fftshift(ifft(fftshift(uddata))));

%creates velocity axis
vfov = 1/(udkv(2)-udkv(1));
venc = vfov/2;
nVE = length(udkv);
dv = vfov/nVE;
vaxis = (-venc):dv:(venc-dv);

figure,
subplot(221),plot(vdkv,abs(vddata),'b-',udkv,abs(uddata),'r-',vdkv,abs(vddata),'b.',udkv,abs(uddata),'r.')
xlabel('k_v (s/cm)'),legend('variable-density','interpolated'),ylabel('Kv-FOCUSING')
subplot(222),plot(vaxis,veldistr,'b-',vaxis,veldistr,'b.')
xlabel('velocity (cm/s)'),

%%%% CONVENTIONAL GRIDDING
[uddata,udkv] = kvfocusing(vddata,vdkv,0,0,1,overgridfactor);
veldistr = abs(fftshift(ifft(fftshift(uddata))));

subplot(223),plot(vdkv,abs(vddata),'b-',udkv,abs(uddata),'r-',vdkv,abs(vddata),'b.',udkv,abs(uddata),'r.')
xlabel('k_v (s/cm)'),legend('variable-density','interpolated'),ylabel('GRIDDING')
subplot(224),plot(vaxis,veldistr,'b-',vaxis,veldistr,'b.')
xlabel('velocity (cm/s)'),