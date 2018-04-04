function test_tet_orientation()
nodeM = -[ 1 0 0; 0 1 0; 0 0 1; 0 0 0];
elem = [1 2 3 4];

N = size(nodeM,1);
SingleTet = elasticTetObj(nodeM, elem);

Y = 100000; % Young's modululs
P = 0.4; % Poisson ratio
rho = 1; % density

SingleTet.SetMaterial( Y, P, rho, 1, 1); % set the first tet to be neohookean
SingleTet.finalize(); % finalize the material of the object

% set the initial state as completely undeformed state with zero energy
Dx = zeros(3*N,1); % displacement field
v = zeros(length(Dx),1);
SingleTet.SetCurrentState(Dx,v);
SingleTet.W


rand('state',1); % Always the same results
fd=@(p) dsphere(p,0,0,0,1);
[p,t]=distmeshnd(fd,@huniform,0.5,1.1*[-1,-1,-1;1,1,1],[]);
nodeM = p;
elem = t;
N = size(nodeM,1);
BallTet = elasticTetObj(nodeM, elem);

Y = 10; % Young's modululs
P = 0.4; % Poisson ratio
rho = 1; % density

BallTet.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tets to linear elasticity
BallTet.finalize(); % finalize the material of the object

% set the initial state as completely undeformed state with zero energy
Dx = zeros(3*N,1); % displacement field
v = zeros(3*N,1); 
BallTet.SetCurrentState(Dx,v);
BallTet.W
end