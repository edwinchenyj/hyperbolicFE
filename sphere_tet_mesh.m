meshname = sprintf('mesh_data%csphere',fs);
load([meshname '.mat'],'node','elem')
elem = elem+1;
p = node;
t = elem;
[p,t] = fixmesh(p,t);
% simpplot(p,t)
% surftri(p,t)
obj = elasticTetObj(p, t);
obj.simple_vis()