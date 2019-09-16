function y = rect(t,tau,xo)
% Generates the rectangular pulse with width tau, center at xo.

% By Davi Marco Lyra Leite (davi@ieee.org) - 02/22/2011

o1 = size(t);
a = t;

y = zeros(o1);
y(abs(a - xo)<(tau/2)) = 1;
y(abs(a - xo)==(tau/2)) = 0.5;
end