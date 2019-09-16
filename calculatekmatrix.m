function [kmatrix,iuf,vuf] = calculatekmatrix(nVE,nphases,viewordering,undersamplingfactor)

if undersamplingfactor ~= 1 && mod(undersamplingfactor/2,1) ~= 0,
    error('undersampling factor must be a power of 2');
end;

%==========================================================================
switch viewordering,

    case {0,1,2,3,11,12,13,21,22,23}, 
        iuf = undersamplingfactor;
        vuf = 1;

    % unfold-style view ordering
    case 4,
        if(undersamplingfactor == 4),
            iuf = 2;
            vuf = 2;
        else
            error('unfold-style view ordering implemented only for undersamplingfactor = 4')
        end;
    otherwise,
        error('Unexpected view-ordering value')
end;

%==========================================================================
kmatrix = nan(nVE,nphases);

for p = 1:nphases,
    for v=1:nVE,
        switch viewordering,

            % [0] same interleaves for all phases (aliasing will be at v = 0 and f =0)
            case 0, 
                kmatrix(v,p) = 1;

            % [1] alternate interleaves between phases (interleaf view-sharing) (aliasing will be at v = 0 and shifted in f)
            case 1,
                kmatrix(v,p) = mod(p-1,undersamplingfactor)+1;
 
            % [2] alternate interleaves between kvs (aliasing will be at v = 0 and f =0)
            case 2,
                kmatrix(v,p) = mod(v-1,undersamplingfactor)+1;

            % [3] alternate interleaves between phases and kvs (aliasing will be at v = 0 and shifted in f)
            case 3,
                kmatrix(v,p) = mod(p+v-2,undersamplingfactor)+1;

            % [4] unfold-style view ordering (aliasing will be shifted in v and f)
            case 4,
                vv=mod(v-1,undersamplingfactor)+1;
                pp=mod(p-1,undersamplingfactor)+1;
                k=abs(vv-pp);
                if k==0,
                    kmatrix(v,p) = 1;
                elseif k==2
                    kmatrix(v,p) = 2;
                end;
                
            case 11,
                ncols = ceil(nphases/undersamplingfactor);
                switch undersamplingfactor,
                    case 1,
                        kmatrixtemp = repmat(1,[nVE ncols]);
                    case 2,
                        kmatrixtemp = repmat([1 2],[nVE ncols]);
                    case 4,
                        kmatrixtemp = repmat([1 3 2 4],[nVE ncols]);
                    otherwise,
                        error('unexpected undersampling factor for view ordering scheme 11.')
                end;
                kmatrix = kmatrixtemp(:,1:nphases);
                
            case 12,
                nrows = ceil(nVE/undersamplingfactor);
                switch undersamplingfactor,
                    case 1,
                        kmatrixtemp = repmat(1,[nrows nphases]);
                    case 2,
                        kmatrixtemp = repmat([1;2],[nrows nphases]);
                    case 4,
                        kmatrixtemp = repmat([1;3;2;4],[nrows nphases]);
                    otherwise,
                        error('unexpected undersampling factor for view ordering scheme 11.')
                end;
                kmatrix = kmatrixtemp(1:nVE,:);
                
            case 13,
                nrows = ceil(nVE/undersamplingfactor);
                ncols = ceil(nphases/undersamplingfactor);
                switch undersamplingfactor,
                    case 1,
                        kmatrixtemp = repmat(1,[nrows ncols]);
                    case 2,
                        kmatrixtemplate = [1 2;
                                           2 1];
                        kmatrixtemp = repmat(kmatrixtemplate,[nrows ncols]);
                    case 4,
                        kmatrixtemplate = [1 3 2 4;
                                           3 2 4 1;
                                           2 4 1 3;
                                           4 1 3 2];
                        kmatrixtemp = repmat(kmatrixtemplate,[nrows ncols]);
                    otherwise,
                        error('unexpected undersampling factor for view ordering scheme 11.')
                end;
                kmatrix = kmatrixtemp(1:nVE,1:nphases);

            %RANDOM VIEW-ORDERING SCHEMES
            case 21,
                krow=floor(1+undersamplingfactor*rand(1,nphases));
                kmatrix = repmat(krow,[nVE,1]);
            case 22,
                kcolumn=floor(1+undersamplingfactor*rand(nVE,1));
                kmatrix = repmat(kcolumn,[1,nphases]);
            case 23,
                kmatrix = floor(1+undersamplingfactor*rand(nVE,nphases));
        end;
    end;
end;