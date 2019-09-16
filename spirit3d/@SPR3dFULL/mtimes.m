function res = mtimes(a,b)
% res = mtimes(FT, x)
%


if a.adjoint  % k -> x ([Nk,Nspr,Nz,Nc] --> [x,y,z,Nc])
    
    res = zeros( a.imSize(1), a.imSize(2), a.imSize(3), a.Nc);
    for cc=1:a.Nc
        
        %         res_cc = zeros( a.imSize(1), a.imSize(2), a.imSize(3));
        %
        kr_z = fftshift( ifft(fftshift( squeeze(b(:,:,:,cc)),3),[],3),3) *sqrt(a.dataSize(3));
        
        for zz=1:a.dataSize(3)
            temp = demodul(kr_z(:,:,zz), a.f0(zz));
            res(:,:,zz,cc) = a.FT{zz}'*temp;
%            res(:,:,zz,cc) = a.FT{zz}'* kr_z(:,:,zz);
            
        end
        
    end
    
    
else % x -> k ([x,y,z,Nc] --> [Nk,Nspr,Nz,Nc]
    res = zeros( a.dataSize(1), a.dataSize(2), a.dataSize(3), a.Nc );
    for cc=1:a.Nc
        
       % im = fftshift( fft(fftshift( squeeze( b(:,:,:,cc)),3),[],3),3) /sqrt(a.dataSize(3));
        
       kr_z = zeros(a.dataSize(1), a.dataSize(2), a.imSize(3));
       for zz=1: a.imSize(3)
           temp = a.FT{zz}* squeeze( b(:,:,zz,cc));
           kr_z(:,:,zz) = demodul( temp, -a.f0(zz));
           %kr_z(:,:,zz) = a.FT{zz}* squeeze( b(:,:,zz,cc));
            
       end
       
       kr_kz = fftshift( fft(fftshift(kr_z,3),[],3),3) / sqrt(a.imSize(3));
        
       for zz=1: a.dataSize(3)
           res(:,a.sampidx(:,zz),zz,cc) = kr_kz(:,a.sampidx(:,zz),zz);
       end
        
        
    end
    
end


  