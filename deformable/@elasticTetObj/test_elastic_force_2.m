function test_elastic_force_2()
clear all
close all


%% ball
rng('default'); 
rng(1) % Always the same results
fd=@(p) dsphere(p,0,0,0,1);
[p,t]=distmeshnd(fd,@huniform,0.4,1.1*[-1,-1,-1;1,1,1],[]);
% pv=[-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
%         1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
%     [p,t]=distmesh2d(@dpoly,@huniform,0.1,[-1,-1; 2,1],pv,pv);
nodeM = 50*p;
elem = t;
% nodeM = [ 0 0 0; 1 0 0; 0 1 0; 0 0 1];
% elem = [1 2 3 4];
N = size(nodeM,1);
M = size(elem,1);
for i = 1:size(elem,1)
    T = [nodeM(elem(i,:),:), ones(4,1)];
    assert(det(T)<0)
end
obj = elasticTetObj(nodeM, elem);

Y = 1000; % Young's modululs
P = 0.4; % Poisson ratio
rho = 1; % density

obj.SetMaterial( Y, P, rho, 1, 0, 0); % set the tets to be neohookean

% obj.SetMaterial( Y, P, rho, 2, 0, 0); % set the tets to be linear

% obj.SetMaterial( Y, P, rho, 3, 0, 0); % set the tets to be stvk

% set the initial state as completely undeformed state with zero energy
Dx = 0e-3*rand(3*N,1); % displacement field
v = zeros(length(Dx),1);
obj.SetCurrentState(Dx);
u0 = [Dx;v]; % the initial state vector
%%

force = obj.ElasticForce;
ep = 1e-3;
rng('shuffle')
normc_fcn = @(m) sqrt(m.^2 ./ sum(m.^2));
dx = normc_fcn(rand(3*N,1));

tic
K = obj.StiffnessMatrix;
toc

obj.SetCurrentState(Dx + ep*dx);
K_new = obj.StiffnessMatrix;
force_new = obj.ElasticForce;

disp('max(df+Kdx)');
disp(max((force_new-force)+ep*K*dx));
disp('max(df+K_newdx)');
disp(max((force_new-force)+ep*K_new*dx));


size(dx)
end