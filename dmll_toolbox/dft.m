% dft.m
%
% neste script são mostrados alguns exemplos da DFT

% André Noll Barreto, 07 jun 2010
N = 64;
Ta = 0.1;
T0 = Ta * N;
fa = 1/T0;
f0 = 1/Ta;

t = ((1:N)-1) * Ta;
t_aux = [((1:N/2)-1), (-N/2 : -1)] * Ta;
f = ([1:N]-1).*fa;

td = 0:N-1;

%
% função porta
T = 1;
g1 = 1.0*(abs(t_aux) <= T);
T=3;
g2 = 0.5*(abs(t_aux) <= T);
G1 = fft(g1);
G2 = fft(g2);

close all

figure(1)
plot(td, g1, td, g2);
title('função rect');
xlabel('n');
ylabel('g_n');

axis([0 N-1 0 2])


figure(2)
plot(td, G1, td, G2);
title('função rect');
xlabel('k');
ylabel('G_k');
axis([0 N-1 -10 10])


%
% função impulso
g1 = 1.0 * (td==0);
% impulso com atraso
g2 = 1.0 * (td==10);
G1 = fft(g1);
G2 = fft(g2);

figure(3)
plot(td, g1, td, g2);
title('função impulso');
xlabel('n');
ylabel('g_n');

axis([0,N-1,0,2]);

figure(4)
plot(td, abs(G1), td, abs(G2));
title('função impulso');
xlabel('k');
ylabel('|G_k|');
axis([0,N-1,0,2]);

figure(5)
plot(td, angle(G1), td, angle(G2));
title('função impulso');
xlabel('k');
ylabel('phi(G_k)');
axis([0,N-1,-pi,pi]);





