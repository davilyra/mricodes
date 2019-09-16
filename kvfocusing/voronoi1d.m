function weight = voronoi1d(ktraj)

%finding w, the pre-weighting function
if(isempty(find(ktraj==0,1))),
    %symmetic sampling
    deltak = ktraj(2:end)-ktraj(1:end-1);
    deltak = [deltak(1);deltak;deltak(end)];
    weight=(deltak(1:end-1)+deltak(2:end))/2;
	weight = weight/min(weight);
    warning('Accurate scaling is not implemented for symmetric kspace coverage!')
else
    %k-space origin is sampled
	ktemp = [ktraj;-ktraj(1)];
	deltak = ktemp(2:end)-ktemp(1:end-1);
	deltak = [deltak(1);deltak];
	weight=(deltak(1:end-1)+deltak(2:end))/2;
	weight = weight/min(weight);
end;