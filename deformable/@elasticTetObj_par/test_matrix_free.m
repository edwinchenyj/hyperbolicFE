function test_matrix_free()
clear all
close all
%% Initialize the object
% a single tet
node = [ 1000 0 0; 0 100 0; 0 0 100; 1 3 5];
elem = [1 2 3 4];
N = size(node,1);
obj = elasticTetObj(node, elem);

Y = 10; % Young's modululs
P = 0.4; % Poisson ratio
rho = 1; % density

obj.SetMaterial( Y, P, rho, 1, 1); % set the first tet to be neohookean
obj.finalize(); % finalize the material of the object


x = 500*rand(3*N,1); % displacement field
v = zeros(3*N,1); 
obj.SetCurrentState(x,v);
fox = obj.ElasticForce();
e1 = []; e2 = []; e3 = [];
tc = 20;
for i = 1:tc
dx = rand(3*N,1);
dx = normc(dx);
ep = 1e-3;

dfodx = obj.ElasticForceDifferential(dx);
Kdx = obj.StiffnessMatrix() * dx;

obj.SetCurrentState(x + ep*dx, v);
foxpdx  = obj.ElasticForce();
directional_derivative = (foxpdx - fox)/ep;
e1 = [e1 max(max(dfodx - directional_derivative))]
e2 = [e2 max(max((Kdx + dfodx)'))]
e3 = [e3 max(max(Kdx + directional_derivative))]

1;
end

plot(1:tc, e1, 1:tc, e2, 1:tc, e3)



end