function rect_horizontal_flip

clear all
close all

fs = filesep;

mesh_shape = 'rect';
maxA = 0.0005;
simulation_type = 'CG';

meshname = sprintf('mesh_data%c%s_maxA_%.d',fs,mesh_shape, maxA);

if exist([meshname '.mat'], 'file') ~= 2
    disp('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end
nodeM = nodeM(:,[2 1]);
elem = elem(:,[1 3 2]);
N = size(nodeM,1);


mesh_shape = 'rect_horizontal';

meshname_new = sprintf('mesh_data%c%s_maxA_%.d',fs,mesh_shape, maxA);

save([meshname_new '.mat'],'nodeM','elem')
        