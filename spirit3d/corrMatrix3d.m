function [AtA,A] = corrMatrix3d(kCalib, kSize)
% (c) Michael Lustig 2007
%
% Modified by Taehoon Shin Mar 2011

[Nx,Ny,Nz,nCoil] = size(kCalib);

A = [];
y = [];

% kSize_xy = [ size(kSize,1) size(kSize,2) ];
% kSize_z = size(kSize,3);

for cc=1:nCoil
    A_cc = [];
    for zz=1:Nz-kSize(3) + 1
        A_zz = [];
        for kk=1: kSize(3)
           tmp = im2col( squeeze(kCalib(:,:,zz+kk-1,cc)), [kSize(1) kSize(2)], 'sliding').';
           A_zz = [ A_zz, tmp];
        end
        %size(tmp)
        A_cc = [ A_cc ; A_zz];
    end
    A = [A, A_cc];
end

%size(A)

AtA = A'*A;
