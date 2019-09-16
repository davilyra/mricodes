function N = noise(s1,s2)
% Calculates the noise amogst two different signals

% By Davi Marco Lyra Leite

s1 = s1(:);
s2 = s2(:);

N = (sum((abs(s1 - s2)).^2));
end