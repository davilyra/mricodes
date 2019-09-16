function E = energyd(fx)
% E = energyd(fx)
% Computes the energy of a signal fx.
% Input:
%    - fx: signal
% Output
%   - E: energy of fx
%
% April 2011
% Written by Davi M. Lyra-Leite <davi@ieee.org>
%
% Modified May 2014

fx =  fx(:);

E = (sum(abs(fx).^2));
end
