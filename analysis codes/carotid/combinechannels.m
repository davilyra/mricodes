function m = combinechannels(x ,dim, pnorm)

x_mag = (sum(abs(x.^pnorm),dim)).^(1/pnorm);
x_sum = sum(x,dim);

phasem = angle(x_sum);

m = x_mag.*exp(1i*phasem);