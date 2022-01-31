fs = filesep;
meshname = sprintf('mesh_data%cbar_example',fs);
load([meshname '.mat'],'nodeM','elem')
% elem = elem;
p = nodeM;
t = elem;
% [p,t] = fixmesh(p,t);
% simpplot(p,t)
% surftri(p,t)
obj = elasticTetObj(nodeM, elem);
obj.simple_vis()