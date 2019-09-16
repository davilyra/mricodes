function sumofsquares = SOS(signal)
% Calculates the Sum of Squares from a given acquired data in parallel
% imaging

% By Davi Marco Lyra Leite (davi@ieee.org) - 09/12/2011


o = size(signal);

sumofsquares = zeros(o(1),o(2));

for i = 1 : o(3)
    sumofsquares = (abs(signal(:,:,i))).^2 + sumofsquares;
end

sumofsquares = sqrt(sumofsquares);

end