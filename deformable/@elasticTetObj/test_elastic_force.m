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
Y = 10000; % Young's modululs
P = 0.1; % Poisson ratio
rho = 1; % density

obj.SetMaterial( Y, P, rho, 2, 0, 0); % set the first tet to be neohookean

normc_fcn = @(m) sqrt(m.^2 ./ sum(m.^2));

x = 0.1 *rand(3*N,1); % displacement field
obj.SetCurrentState(x);
figure
tetramesh(elem, obj.node)
force = obj.ElasticForce;

ep = 1e-11;
dx = normc_fcn(rand(3*N,1));

K = obj.StiffnessMatrix;
obj.SetCurrentState(ep*dx + x);

force_new = obj.ElasticForce;
disp('f directional derivative')
disp((force_new-force)/ep)
disp('-Kdx')
disp(-K*dx)

disp('max(f+Kdx)')
disp(max((force_new-force)/ep+K*dx))
end