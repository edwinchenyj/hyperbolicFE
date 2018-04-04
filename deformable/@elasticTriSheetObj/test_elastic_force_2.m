function test_elastic_force_2()
clear all
close all


%% ball
rng('default'); 
rng(1) % Always the same results
% nodeM = [1 0; 0 1; 0 0];
% elem = [1 2 3];
% N = size(nodeM,1)
% M = size(elem,1)


    nodeM = [                % list of xy "node" coordinates
    1, 0
    1, 1
    0, 1
    0, 0                
        ] ;
    
    edge = [                % list of "edges" between nodes
        1, 2                % outer square 
        2, 3
        3, 4
        4, 1
        ] ;

%    [nodeM,etri, ...
%     elem,tnum] = refine2(nodeM,edge) ;
   
hfun = +.5 ;     
      [nodeM,etri, ...
    elem,tnum] = refine2(nodeM,edge,[],[],hfun) ;
trimesh(elem,nodeM(:,1),nodeM(:,2));
N = size(nodeM,1);
M = size(elem,1);
% for i = 1:size(elem,1)
%     T = [nodeM(elem(i,:),:), ones(3,1)];
%     assert(det(T)<0)
% end
obj = elasticTriObj(nodeM, elem);

Y = 100; % Young's modululs
P = 0.45; % Poisson ratio
rho = 1; % density

obj.SetMaterial( Y, P, rho, 1, 0, 0); % set the tets to be neohookean

% set the initial state as completely undeformed state with zero energy
Dx = 0*rand(2*N,1); % displacement field
obj.SetCurrentState(Dx);
% u0 = [Dx]; % the initial state vector
%%

force = obj.ElasticForce;
ep = 1e-8;
rng('shuffle')
dx = normc(rand(2*N,1));
% 
% 
% disp('ElasticForceDifferential')
% tic
% df = obj.ElasticForceDifferential(dx);
% toc

tic
K = obj.StiffnessMatrix;
toc

% obj.SetCurrentState(Dx + ep*dx);
obj.SetCurrentState(ep*dx);
K_new = obj.StiffnessMatrix;
force_new = obj.ElasticForce;

disp('max(df+Kdx)');
disp(max((force_new-force)/ep+K*dx));
% disp(max((force_new-force)/ep-df));
% disp(max(df+K*dx))

size(dx)
end