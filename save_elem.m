function save_elem
%% set parameters, flags, and meshes
close all

mesh_shape = 'four_elem';
fs = filesep;

elem = [1 2 3; 2 4 3; 2 5 4; 5 6 4];
nodeM = [0 0; 1 0; 0 1; 1 1; 2 0; 2 1]/2;



meshname = sprintf('mesh_data%c%s',fs,mesh_shape);

save([meshname '.mat'],'nodeM','elem')