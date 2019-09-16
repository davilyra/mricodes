function [dr,usercv,hdr] = rawloadHD_jfn(fn,frsize,ep,ishosp, rhbline, coils, slices, echoes)
% function [dr,usercv,hdr] = rawloadHD_jfn(fn,frsize,ep,ishosp, rhbline, coils, slices, echoes)
%
% Similar to rawloadHD, except returns multi-dimensional raw data array.
% In addition, picks out data for desired coils, slices, and echoes.
% This is useful for very large Pfiles on computers with limited RAM.
% 
% Please check to make sure this function does what you want.
%
%
% See also:  pfileinfo_jfn.m.
%
% INPUTS:
%   fn, frsize, ep, ishosp -- se rawloadHD.m
%   rhbline - 0 or 1
%   coils   - return data array containing data for these coils only.
%   slices  - return data array containing data for these slices only.
%   echoes  - return data array containing data for these echoes only.
%
% OUTPUT:
%   dr         - array containing complex raw data
%                Dimensions are frsize x rhnframes x ncoils x nslices x nechoes.
%                For example, d(:,:,2,1,3) is the 2D raw data for coil 2,
%                slice 1, echo 3.
%   usercv     - rhuser1-19
%   hdr        - ?
%
%
% jfn, 2006
% $Id: rawloadHD_jfn.m,v 1.37 2007-01-23 23:41:02 jfnielse Exp $


if (nargin < 4) | (length(ishosp)==0)
	ishosp=0;
end;

%
%  **** For more information, the P-file header is well documented in 
%	the LX ESE users manual.
%

if (nargin < 1)	% Get last p-file if function exists.
	if (exist('lastpfile'))
	    fn = lastpfile;
	    tt = sprintf('No file given... using most recent file, %s.',fn);
	    disp(tt);
	else
		error('You must specify a P-file');
	end;
end;

% default to normal precision
if (nargin < 3) | (length(ep)==0)
  ep = 0;
end

fip = fopen(fn,'r'); % Took away 'b' with Excite upgrade...
if fip == -1
  tt = sprintf('File %s not found\n',fn);
  dispsafe(tt);
  return;
end

ver = fread(fip,1,'float');
if 0 % jfn
dispsafe('Raw file version:');
end


tothdrsize = 66072;  % 60464 in EX  %39940 in LX;
hdrrecsize = 2048;
daqtabsize = 20480;
if ((ver > 7) & (ishosp==0))
	daqtabsize = 2*20480;	% Fixed?!
end;
othhdrsize = 4096+4096+daqtabsize+2052+2052+2048+1000;  % 1000 is empirical
examhdrsize = 1036;	% Was 1040, changed Oct 27/03 
serieshdrsize = 984; %1028;  Could this be the 44-byte difference??
imagehdrsize = 1044;

% These should add to tothdrsize.
hdrrecsize+othhdrsize+examhdrsize+serieshdrsize+imagehdrsize;

% ===========================================================
% Extract Date and Time from header.
% ===========================================================

nstr = 40;
%nstr = 44;

hstr = fread(fip,nstr,'char');
hstr = char(hstr(:))';
if 0 % jfn
dispsafe(hstr);
end;

dstr = sprintf('%c',hstr(13:22));
if (dstr(7:8)=='10')
	dstr = [dstr(1:6) '200' dstr(9)];
end;
tstr = sprintf('%c',hstr(23:30));

if 0 % jfn
disp(' ');
disp('============================================================');
disp('	RawloadHD --> HD Version.				');
disp('============================================================');
dt = sprintf('Raw File:  %s ',fn);
dispsafe(dt);
dt = sprintf('Data Acquisition Date/Time: %s, %s',dstr,tstr);
dispsafe(dt);
end


% ===========================================================
% Read part of the header up to and including rhuser0
% ===========================================================
hdr = fread(fip,88,'int16');

if (nargin < 2) | (length(frsize)==0)
  frsize = hdr(19);
  dt = sprintf('Detected Framesize of %d',frsize);
else
  dt = sprintf('Specified Framesize of %d',frsize);
end
if 0 % jfn
dispsafe(dt);
end
if (nargin < 4) | (length(ishosp)==0)
	ishosp=0;
end;

rhnframes=hdr(16);
nslices=hdr(13);
nechoes=hdr(14);
nex=hdr(15);

dt = sprintf('%d frames, %d slice(s), %d echo(es), %d NEX',rhnframes,nslices,nechoes,nex);
if 0 % jfn
dispsafe(dt);
end

rawsize = 32768-sign(hdr(37))*32768 + hdr(37)+65536*hdr(38);
	% rawsize is a 32-bit value, but stored as 2 16-bit words.
	% The above is the conversion.
bytes_per_sample = 4;
nframes = floor(rawsize / frsize / bytes_per_sample +0.5);
dt = sprintf('%d Total Frames expected',nframes);
if 0 % jfn
dispsafe(dt);
end

clear hdr;


% ===========================================================
% Read rhuser1-rhuser19
% ===========================================================
usercvsize = 19*4;
usercv = fread(fip,usercvsize/4,'float');


% ===========================================================
% Read rest of RDB header.
% ===========================================================
hdr = fread(fip,(hdrrecsize-nstr-usercvsize-176)/4,'int32');

prescanpars = hdr(30:37);	
tt = sprintf('Prescan Parameters:   AX=%ld  TG=%d  R1=%d  R2=%d',prescanpars(8),prescanpars(7),prescanpars(5),prescanpars(6));
if 0 % jfn
dispsafe(tt);
end

fread(fip,11,'float');

% ===========================================================
% Read Next Few headers.
% ===========================================================
hdr = fread(fip,othhdrsize+examhdrsize+serieshdrsize,'char');


% ===========================================================
% Read Image Header, and interpret/display some Info.
% ===========================================================

hdr = fread(fip,floor((imagehdrsize)),'uint8');

% ======= TR, TE =========
tr = bytes2int32(hdr(201:204));
te = bytes2int32(hdr(209:212));
tt = sprintf('TR=%8.3f ms,  TE=%8.3f ms',tr/1000,te/1000);
%dispsafe(tt);	 TR/TE NOT WORKING IN EXCITE HEADER YET!!

% ====== PSD Name ========
tt = sprintf('PSD = %s',char(hdr(321:353))');
%dispsafe(tt);	PSD NOT WORKING WITH EXCITE HEADER YET!!


% ====== Trigger Delay =======
tdel = bytes2int32(hdr(233:236));
hr = bytes2int16(hdr(231:233));
tt = sprintf('Cardiac Info:  Rate = %8.3f bpm  Delay = %8.3f ms',hr,tdel/1000);
%dispsafe(tt);	CARDIAC INFO NOT WORKING WITH EXCITE HEADER YET!!
if 0 % jfn
dispsafe('============================================================');
end


% ===========================================================
% Read all the data
% ===========================================================

% disabled -- jfn
if 0

fseek(fip,tothdrsize,-1);

if ep == 1 
  dr = fread(fip,nframes*frsize*2,'int32');
  dd = fread(fip,inf,'int32');
else
  dr = fread(fip,nframes*frsize*2,'int16');
  dd = fread(fip,inf,'int16');
end;
ll = length(dr);
dr = dr(1:2:ll)+i*dr(2:2:ll);
nframes = floor(ll/2/frsize);
ld = length(dd);
tt = sprintf('%d complex samples discarded at end of file',ld/2);
disp(tt);
dr = dr(1:nframes*frsize);
dr = reshape(dr,frsize,nframes);

end

%% BEGIN additions by jfn

% reshape dr

%rhbline = 1;

sliceres = (frsize+rhbline);
coilres = sliceres*nslices;

%d = reshape(dr, frsize, frsize+rhbline, nslices, 4);

%for coil = 1 : 4
%	for slice = 1 : nslices
%		offset = (coil-1)*coilres + (slice-1)*sliceres + 1;
%		d(:,:,coil,slice) = dr(:,offset:(offset+frsize-1));
%	end
%end



if 0
% print some useful info to screen
fprintf(1, '\nrawloadHD_jfn: \n');
fprintf(1, '\tfrsize    %d\n', frsize);
fprintf(1, '\trhnframes %d\n', rhnframes);
for cvcnt = 1 : 11
	fprintf(1, '\trhuser%d = %f\n', cvcnt, usercv(cvcnt));
end;
end;


if 1
% calculate sizes of data chunks
echores = frsize  * (rhnframes+rhbline);   % number of data points per echo
sliceres = nechoes * echores;              % number of data points per slice
coilres  = nslices * sliceres;             % number of data points per receiver coil

% read data from file into array
% read only the desired coils,  slices, and echoes
dr = zeros(frsize, rhnframes, length(coils), length(slices), length(echoes));
for coilind = 1 : length(coils)
	coil = coils(coilind);
	for sliceind = 1 : length(slices)
		slice = slices(sliceind);
		for echoind = 1 : length(echoes)
			echo = echoes(echoind);
			%offsetres = (coil-1)*coilres+ (slice-1)*sliceres + (echo-1)*echores + frsize*rhbline + rhnframes/2 - frsize+rhnframes;
			offsetres = (coil-1)*coilres+ (slice-1)*sliceres + (echo-1)*echores + frsize*rhbline ;
			offsetbytes = 4 * offsetres ;    % 2 short ints per data sample
			fseek(fip, tothdrsize + offsetbytes, -1);
			beg = ftell(fip);
			d = fread(fip, 2*frsize*rhnframes, 'int16');
			d = complex(d(1:2:end), d(2:2:end));
			%size(d)
			d = reshape(d, frsize, rhnframes);
			dr(:, :, coilind, sliceind, echoind) = d;
			stop = ftell(fip);
			%stop-beg
		end;
	end;
end;
end

fclose(fip);

% END of jfn's additions

return;










function intout = bytes2int32(bytes)
%
%	Convert 4 bytes to an int32
%
intout = (((bytes(1)*256)+bytes(2)*256)+bytes(3)*256)+bytes(4);
if (intout >= 32768*65536)
	intout = intout-65536*65536;
end;



function intout = bytes2int16(bytes)
%
%	Convert 4 bytes to an int32
%
intout = (bytes(1)*256)+bytes(2);
if (intout >= 32768)
	intout = intout-65536;
end;





function notused

d = zeros(frsize, rhnframes, ncoils, nslices, nechoes);
for coil = 1 : ncoils
	% temporary fix for echo mis-centering along ky
	%deltaoffset = 0;
	for slice = 1 : nslices
		for echo = 1 : nechoes
			%offset = (coil-1)*coilres + (slice-1)*sliceres + (echo-1)*echores + rhbline*frsize + 1; 
			offset = (coil-1)*coilres/frsize + (slice-1)*sliceres/frsize + (echo-1)*echores/frsize + rhbline; 
			%dtmp = dr(offset : (offset+echores-1));
			%d(:,:,coil,slice,echo) = reshape(dtmp,frsize,rhnframes);
			d(:,:,coil,slice,echo) = dr(:,(offset+1):(offset+1+rhnframes-1));
			%viewoffset = ((coil-1)*coilres + (slice-1)*sliceres + (echo-1)*echores)/frsize + 1 + deltaoffset; 

			% temporary fix for echo mis-centering along ky
			if 0
				if viewoffset + rhnframes -1 > nframes
					viewoffset = nframes - rhnframes + 1;
				end;
				d(:,:,coil,slice,echo) = dr(:,viewoffset:(viewoffset+rhnframes-1));
			end

			%deltaoffset = deltaoffset + 1;
		end;
	end;
end;

