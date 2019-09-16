function pcref = calcpcref(kv,kvdata)

% finds central kv-samples
[~,idx] = sort(abs(kv));
kv0idx = idx(1);
kv0 = kv(kv0idx);
kv1 = kv(kv0idx+1);
kvdata0 = kvdata(kv0idx);
kvdata1 = kvdata(kv0idx+1);

%calculates phase contrast velocity
vfov = 1/(kv1-kv0);
pcref = angle(kvdata0.*conj(kvdata1));
if(pcref==pi), pcref = 0; end;
pcref = pcref/pi * vfov/2;