function test_elastic_force()
clear all
close all

node = rand(4,3); % set random undeform nodal positions
node = [ 1 0 0; 0 1 0; 0 0 1; 0 0 0];
elem = [1 2 3 4];
N = size(node,1);
G = elasticTetObj.G;
Dm = node' * G;

% need to make sure the random nodal positions generated creates a tet with
% positive volume
if det(Dm) < 0
    node(4,:) = node(4,:) + 2 * Dm(:,1)'; % move node #4 to the other side along the vector from node #4 to node #1 (too lazy to do reflection...)
end
Dm = node' * G;
assert(det(Dm) > 0)

tetramesh(elem, node)

obj = elasticTetObj(node, elem);

vol = obj.W;
% assert(vol > 0)
Y = 1; % Young's modululs
P = 0.1; % Poisson ratio
rho = 1; % density

obj.SetMaterial( Y, P, rho, 1, 1); % set the first tet to be neohookean

% obj.mu = 0;
obj.finalize(); % finalize the material of the object


x = 0.0 *rand(3*N,1); % displacement field
% x(1) = x(1) + 0.5;
% x(4) = x(4) + 0.5;
% x(7) = x(7) + 0.5;
% x(10) = x(10) + 0.5;
v = zeros(3*N,1);
obj.SetCurrentState(x,v);
figure
tetramesh(elem, obj.node)
Ds = obj.Ds;
Dm = obj.Dm;
DmINV = obj.DmINV;
Ds - Dm
Ds * inv(Dm)
disp('F')
disp(obj.F)
disp('P')
P = obj.Stress(1);
disp(P)
disp('force')
force = obj.ElasticForce;
disp(force);

ep = 1e-3;
dx = normc(rand(3*N,1));

dxNode = reshape(dx,3,obj.N)'
t = 1;
dDs = dxNode(obj.elem(t,:),:)' * G;

dF = dDs * obj.DmINV(3*(t-1)+1:3*t,:)
%   dF =  dxNode' * G * DmINV


disp('Stress Differential')
disp(obj.StressDifferential(1,dF))
disp('ElasticForceDifferential')
disp(obj.ElasticForceDifferential(dx))
K = obj.StiffnessMatrix;
obj.SetCurrentState(ep*dx, v);
% disp('new F')
% disp(obj.F)
% disp('new P')
P_new = obj.Stress(1);
% disp(P_new)
disp('directional_derivative')
disp((P_new - P)/ep)

force_new = obj.ElasticForce;
disp('f directional derivative')
disp((force_new-force)/ep)
disp('-Kdx')
disp(-K*dx)

disp('max(f+Kdx)')
disp(max((force_new-force)/ep+K*dx))
end