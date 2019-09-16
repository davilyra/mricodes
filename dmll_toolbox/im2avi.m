function varargout=im2avi(ext, imdir, scale, framerate, filename, playflag)
% im2avi converts image sequenses to avi video
%
% SYNTAX
%
% Inputs: imdir: image sequence directory
% ext: the extension name of the image such as 'jpg', 'tif',
% scale: image resize, like [320 400] or 0.9
% framerate: avi video frame rate
% filename: save avi as
% playflag: play the avi video; if it is 0, no play, if >0,
% play the avi 'playflag' times.
%
% EXAMPLE: im2avi(ext, imdir, scale, framerate, filename, playflag)
%
% NOTES: based on im2avi (author: Zhe Wu @ Univ of Rochester)
% Wenbin, 09-May-2007

warning('im2avi:warning','Do NOT close or open any image window during the process!\n')

filearray=dir([imdir filesep '*.' ext]);
s=size(filearray,1);

frameind=0;
mv =struct('cdata',{}, 'colormap', {});
figure, h =gcf;

for i=1:s
frameind=frameind+1;
imgname=[imdir filesep filearray(i).name];
im=imread(imgname) ;
im=imresize(im, scale);
imshow(im);
mv(frameind)=getframe(h);
end
close(h)

%movie2avi(mv, [imdir filesep filename '.avi'], 'fps', framerate);
movie2avi(mv, [imdir filesep filename '.avi'], 'fps', framerate, 'compression', 'None');

if nargout >0
varargout{1} =mv;
end

if playflag>0
movie(mv, playflag);
end