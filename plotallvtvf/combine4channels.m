function m = combine4channels(m1,m2,m3,m4)
% usage: m = combine4channels(m1,m2,m3,m4);
%
% Combines up to 4 channels. If using less than 4 channels,
% simply set the extra m's to 0.
%
% example: m = combine4channels(m1,m2,m3,0);
im1_mag = abs(m1); im1_phase = angle(m1);
im2_mag = abs(m2); im2_phase = angle(m2);
im3_mag = abs(m3); im3_phase = angle(m3);
im4_mag = abs(m4); im4_phase = angle(m4);

% % im#_mag is magnitude of coil #, im#_phase is phase of coil #
% % combine mag
absm = sqrt( im1_mag.^2 + im2_mag.^2 + im3_mag.^2 + im4_mag.^2);
% % combine phase
im_phase_real_comb = (im1_mag.^2) .* cos(im1_phase) + ...
                     (im2_mag.^2) .* cos(im2_phase) + ...
                     (im3_mag.^2) .* cos(im3_phase) + ...
                     (im4_mag.^2) .* cos(im4_phase);
im_phase_imag_comb = (im1_mag.^2) .* sin(im1_phase) + ...
                     (im2_mag.^2) .* sin(im2_phase) + ...
                     (im3_mag.^2) .* sin(im3_phase) + ...
                     (im4_mag.^2) .* sin(im4_phase);     
phasem = angle(im_phase_real_comb+1i*im_phase_imag_comb);
m = absm.*exp(1i*phasem);