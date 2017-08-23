% generate the path of the main folder that contains all the
% subfolders/object classes
foldername = fileparts(which(mfilename));
cd(foldername)

% generate the pathes to all the subfolders
p = genpath(foldername);

% add everything to the search path
addpath(p);

% remove old sim
% rmpath([foldername '/old sim']);