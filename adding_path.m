% generate the path of the main folder that contains all the
% subfolders/object classes
foldername = fileparts(which(mfilename));
cd(foldername)

% generate the pathes to all the subfolders
p = genpath(foldername);

% add everything to the search path
addpath(p);

% remove sim
sim_p = genpath([foldername filesep 'sim_data']);
% rmpath([foldername filesep 'sim_data']);
rmpath(sim_p);
rmpath([foldername filesep '.git']);
rmpath([foldername filesep 'codegen']);
rmpath([foldername filesep 'mex_files' filesep 'codegen']);
rmpath([foldername filesep 'mex_files' filesep 'codegen' filesep 'mex']);
