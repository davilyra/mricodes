function V = ser(signal1,signal2)
% Calculates the Signal-Noise Ratio in dB for a given data

% By Davi Marco Lyra Leite (davi@ieee.org)

V = 10*log10(sum((abs(signal1(:))).^2)/(sum((abs(signal1(:)-signal2(:))).^2)));
end