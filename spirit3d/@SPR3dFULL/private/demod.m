function raw_demod = demod( raw, df )
% raw_demod = demod( raw, df)
%
% Input
%   raw - [Nk,Nspr] spiral rawdata
%   df - Off-resonance [Hz]
%
% Output
%   raw_demod - [ Nk,Nspr] demodulated rawdata by df
%
%      Taehoon Shin Mar 2009

Nk = size(raw,1);
Nspr = size(raw,2);

tmap = 4e-6*[ 0: Nk-1];
phmap = exp( -2*j*pi*df*tmap );
raw_demod = zeros(Nk,Nspr);
for ss=1:Nspr
    raw_demod(:,ss) = reshape(raw(:,ss),1,Nk) .* phmap;
end

% figure; plot( 1:Nk, abs( raw(:,1,6)), 1:Nk, abs( raw_demod(:,1,6))); 
% 
% figure; plot( 1:Nk, angle( raw(:,1,6)), 1:Nk, angle( raw_demod(:,1,6))); 
