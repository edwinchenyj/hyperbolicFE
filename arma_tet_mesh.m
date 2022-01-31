% save tetgen mesh from read_tetgenmesh
name = 'arma_6';
fs = filesep;
meshname = sprintf('mesh_data%c%s',fs,name);
% if exist([meshname '.mat'], 'file') ~= 2
    [p,t] = read_tetgenmesh(meshname);
    t = t + 1;
    t = t';
    p  = p(1:3,:)';
% elem = elem;
% p = node;
% t = elem;
% [p,t] = fixmesh(p,t);
% simpplot(p,t)
% surftri(p,t)

obj = elasticTetObj(p, t);

max_z = max(p(:,3));
free_ind = reshape(repmat(p(:,3) < max_z - 0.12,1,3)',[],1);

Y = 1e5;
P = 0.45;
rho = 1000;
a = 0;
b = 0;
obj.SetMaterial( Y, P, rho, 1, a, b);
Dx = 0*rand(obj.Dim*obj.N,1);
obj.SetCurrentState(Dx);

v = zeros(length(Dx),1);
u = [Dx; v];

obj.gravity_on = true;
obj.calculateGravity;

obj.simple_vis()
dt = 5e-3;
% v = VideoWriter('arma_sim');
% open(v)
for i = 1:200

obj.simple_vis()
if i == 51 
    1;
end
u = SemiImplicit(dt, u, obj,~free_ind);

% writeVideo(v,getframe(gcf))
end

% close(v)