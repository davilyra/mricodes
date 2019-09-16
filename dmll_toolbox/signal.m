function S = signal(s1)
% Calculates the power of a signal

% By Davi Marco Lyra Leite

s1 = s1(:);

S = (sum((abs(s1)).^2));
end