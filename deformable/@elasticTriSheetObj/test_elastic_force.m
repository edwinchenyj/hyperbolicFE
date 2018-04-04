function test_elastic_force()
clear all
close all

node = rand(3,2); % set random undeform nodal positions
node = [ 1 0; 0 1; 0 0];
elem = [1 2 3];
N = size(node,1);
G = staticTriMesh.G;
Dm = node' * G;

% need to make sure the random nodal positions generated creates a tet with
% positive volume
if det(Dm) < 0
    node(3,:) = node(3,:) + 2 * Dm(:,1)'; % move node #4 to the other side along the vector from node #4 to node #1 (too lazy to do reflection...)
end
Dm = node' * G;
assert(det(Dm) > 0)

% tetramesh(elem, node)
trimesh(elem,node(:,1),node(:,2))

obj = staticTriMesh(node, elem);

vol = obj.W;
% assert(vol > 0)
Y = 10000000; % Young's modululs
P = 0.1; % Poisson ratio
rho = 1; % density

obj.SetMaterial( Y, P, rho, 1, 1); % set the first tet to be neohookean

x = 0.0 *rand(2*N,1); % displacement field
% x(1) = x(1) + 0.5;
% x(4) = x(4) + 0.5;
% x(7) = x(7) + 0.5;
% x(10) = x(10) + 0.5;
obj.SetCurrentState(x);
figure
Ds = obj.Ds;
Dm = obj.Dm;
DmINV = obj.DmINV;
Ds - Dm;
Ds * inv(Dm);
disp('F')
disp(obj.F)
disp('P')
P = obj.Stress(1);
disp(P)
disp('force')
force = obj.ElasticForce;
disp(force);
ep_list = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
for i = 1:length(ep_list)
ep = ep_list(i);
dx = normc(rand(2*N,1));

dxNode = reshape(dx,2,obj.N)';
t = 1;
dDs = dxNode(obj.elem(t,:),:)' * G;

dF = dDs * obj.DmINV(2*(t-1)+1:2*t,:);
%   dF =  dxNode' * G * DmINV
% 
% 
% disp('Stress Differential')
% disp(obj.StressDifferential(1,dF))
% disp('ElasticForceDifferential')
% disp(obj.ElasticForceDifferential(dx))
K = obj.StiffnessMatrix;
obj.SetCurrentState(ep*dx);
% disp('new F')
% disp(obj.F)
% disp('new P')
P_new = obj.Stress(1);
% disp(P_new)
% disp('directional_derivative')
% disp((P_new - P)/ep)

force_new = obj.ElasticForce;
% disp('f directional derivative')
% disp((force_new-force)/ep)
% disp('-Kdx')
% disp(-K*dx)

disp('max(f+Kdx)')
disp(max((force_new-force)/ep+K*dx))
end
end