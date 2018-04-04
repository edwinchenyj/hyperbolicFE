function test_object_deep_copy()

%%
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



obj = elasticTetObj(node, elem);

vol = obj.W;
assert(vol > 0)
Y = 1; % Young's modululs
P = 0.1; % Poisson ratio
rho = 1; % density

obj.SetMaterial( Y, P, rho, 1, 1); % set the first tet to be neohookean

% obj.mu = 0;
obj.finalize(); % finalize the material of the object



obj2 = elasticTetObj(obj);


x1 = 0.5 *rand(3*N,1); % displacement field
v = zeros(3*N,1);
obj.SetCurrentState(x1,v);

x2 = 0.5 *rand(3*N,1); % displacement field
v = zeros(3*N,1);
obj2.SetCurrentState(x2,v);

%%
tetNodes1 = obj.GetX; tetNodes1 = reshape(tetNodes1,3,N)'; % re-organize the nodes into (N by 3)


figure('Name','tet1')
tetramesh(elem,tetNodes1);

tetNodes2 = obj2.GetX; tetNodes2 = reshape(tetNodes2,3,N)'; % re-organize the nodes into (N by 3)

figure('Name','tet2')
tetramesh(elem,tetNodes2);

end