function avb = amsv (ama,svx)

% This function was taken from the 'library_vectorization' folder of the
% MATLAB File Exchange submission
%
%   http://www.mathworks.com/matlabcentral/fileexchange/27826
%
% by Talal Rahman and Jan Valdman.


% ama: ama(1:nx,1:ny,1:nz)
% svx: svx(1:ny,1)
% avb: avb(1:nx,1,1:nz)

[nx,ny,nz] = size(ama);

avx = svx(:).';
avx = avx(ones(nx,1),:,ones(nz,1));

avb = ama .* avx;
avb = sum(avb,2);

return