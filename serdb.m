function x = serdb(signal,estimate)

signalenergy= sum(abs(signal(:)).^2);
errorenergy= sum(abs(signal(:)-estimate(:)).^2);

x = 10*log10(signalenergy/errorenergy);