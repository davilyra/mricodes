function mwdea = grid2_deap(d,k,n)
% function m = grid2(d,k,n)
% Performs the gridding and the deapodization process
%
% Inputs
%   - d: k-space data
%   - k: k-trajectory, scaled -1 to 1 (2X gridding)
%   - n: image size
%
% Outputs:
%   - mwdea: deapodized image in the image-space
%
% Original code of the gridding by Krishna Nayak, PhD,
% University of Southen California, December 1st 2009
%
% This version by Davi Marco Lyra Leite (davii@ieee.org), Apr 23 2013

% convert to single column
d = d(:);
k = k(:);

% convert k-space samples to matrix indices
nx = ((n/2+1) + n*real(k));
ny = ((n/2+1) + n*imag(k));

% zero out output array
m = zeros(2*n,2*n);

% loop over samples in kernel
for lx = -2:2,
	for ly = -2:2,
		% find nearest samples
		nxt = round(nx + lx);
		nyt = round(ny + ly);
        nxt2 = round(nx + lx) + 0.5;
        nyt2 = round(ny + ly) + 0.5;

		% compute weighting for triangular kernel
		kwx = max(1-abs(nx-nxt),0);
		kwy = max(1-abs(ny-nyt),0);
        kwx2 = max(1-abs(nx-nxt2),0);
        kwy2 = max(1-abs(ny-nyt2),0);

		% map samples outside the matrix to the edges
		nxt = max(nxt,2); nxt = min(nxt,n);
		nyt = max(nyt,2); nyt = min(nyt,n);
        nxt2 = max(nxt2,1); nxt2 = min(nxt2,n);
        nyt2 = max(nyt2,1); nyt2 = min(nyt2,n);

		% use sparse matrix to turn k-space trajectory into 2D matrix
		m = m + sparse([2*nxt2-1 2*nxt-1 2*nxt2-1 2*nxt-1],[2*nyt2-1 2*nyt-1 2*nyt-1 2*nyt2-1],[(d.*kwx2.*kwy2) (d.*kwx.*kwy) (d.*kwx2.*kwy) (d.*kwx.*kwy2)],n*2,n*2);
	end;
end;

% zero out edge samples, since these may be due to samples outside
% the matrix
m(:,1) = 0; m(:,2*n) = 0;
m(1,:) = 0; m(2*n,:) = 0;


%% Dapodization:
x = -1:(1/n):(1-1/n);
y = x;
[xx,yy] = meshgrid(x,y);

finvx = pi.*(sinc(xx/1.25)).^2;
finvy = pi.*(sinc(yy/1.25)).^2;

finv = finvx.*finvy;

% Second element weighted data
mwdea = ift(m)./finv;

% delivering the cropped data
mwdea = mwdea((round(n/2 + 1)):(round(3*n/2)),(round(n/2 + 1)):(round(3*n/2)));
