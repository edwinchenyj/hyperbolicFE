function test_elastic_force_2()
clear all
close all


D = 10;
[x,y] = meshgrid(1:D,1:D);
tri = delaunay(x,y);
z = peaks(D)/10;

nodeM = [x(:), y(:), zeros(size(x(:)))];
nodeM = [x(:), y(:), z(:)];

Adj = sparse(size(nodeM,1),size(nodeM,1));
elem = [];
for t = 1:size(tri,1)
    node1 = tri(t,1);
    node2 = tri(t,2);
    node3 = tri(t,3);
    if ~Adj(node1,node2)
        Adj(node1,node2) = 1;
        Adj(node2,node1) = 1;
        elem = [elem; node1 node2];
    end
    if ~Adj(node2,node3)
        Adj(node2,node3) = 1;
        Adj(node3,node2) = 1;
        elem = [elem; node2 node3];
    end

    if ~Adj(node1,node3)
        Adj(node3,node1) = 1;
        Adj(node1,node3) = 1;
        elem = [elem; node3 node1];
    end

end
% nodeM = [0 0 0; 1 1 1];
% elem = [1 2];
% initialize the object
% a line of strand elements
% initialize n particles on a vertial chain with seperation l(i)
n = size(nodeM,1);
% fixing a random seed
% initialize m strands
m = size(elem,1);
% stiffness
k = 10000*rand(m,1);
% linear damping constants
k_damping = 0 * rand(m,1);
% each row contains the indices of the particles connected by that spring

% initial positions
q = nodeM';
q = q(:);

% constructor 
% obj = springs_particles(q,v,springs,k,k_damping,l,mass)
obj = elasticStrandObj(nodeM,elem);
rho = 1;
for s = 1:m
obj.SetMaterial(k(s),rho,s,1);
end
obj.finalize();


% set the initial state as completely undeformed state with zero energy
Dx = zeros(3*n,1);
v = zeros(length(Dx),1);
obj.SetCurrentState(Dx,v);
Dx = Dx(4:end); v = v(4:end); % fix the first node
u0 = [Dx;v];

%%

force = obj.ElasticForce;
ep = 1e-5;
rng('default')
rng(1);
dx = normc(rand(3*n,1));


disp('ElasticForceDifferential')
tic
df = obj.ElasticForceDifferential(dx);
toc

tic
K = obj.StiffnessMatrix;
toc

obj.SetCurrentState(ep*dx, v);

force_new = obj.ElasticForce;

disp('max(df+Kdx)');
disp(max((force_new-force)/ep+K*dx));
disp(max((force_new-force)/ep-df));
disp(max(df+K*dx))

size(dx)
end