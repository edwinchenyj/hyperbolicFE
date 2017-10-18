function test_DG_b

clear all
close all

draw = true;
rerun_flag = true;
save_state = true;
draw = true;

dt = 1/120;
tsteps = 120*10;

fs = filesep;

mesh_shape = 'triangle';

constraints = 1; % types of constraint
% 1: free

deformation_mode_number = 1;

Y = 10; % Young's modululs
P = 0.48; % Poisson ratio
rho = 1; % density
a = 0.0; % rayleigh damping
b = 0.00;

elem = [1 2 3; 2 4 3];
nodeM = [-1 -1; 1 0; 0 1; 1 1]/2;
N = 4;

% construct triangular mesh object
obj = elasticDGTriObj(nodeM, elem);


obj.SetMaterial( Y, P, rho, 1:size(elem,1), 2); % set the tri to linear
Dx = 0*rand(2*N,1); % displacement field. set zero for rest state
obj.SetCurrentState(Dx);
%
M = obj.M;
K = obj.StiffnessMatrix;
%         K = K(~ind_fix,~ind_fix); % extract the non-fixed part

[V,D] = eig(full(K),full(M));
[low_eig, permutation_indices] = sort(diag(D));
V = -V(:,permutation_indices);
firstMode = V(:,4)/2;

DGN = size(obj.DGnodeM,1);

TR = triangulation(obj.DGelem, obj.DGnodeM);
E = edges(TR);

%%
% deformation for the initial condition
deformation_scale_factor = 10;
deformation_mode = V(:,3 + deformation_mode_number)/deformation_scale_factor;

node = obj.node;
positionsM = node';
positionsM = positionsM(:);
positions = positionsM;
externalGravity  = zeros(size(positions));
nFixed = 0;
indLogical = true(size(positions));

DGnode = obj.DGnode;
DGpositionsM = DGnode';
DGpositionsM = DGpositionsM(:);
DGpositions = DGpositionsM;
DGexternalGravity  = zeros(size(DGpositions));
DGindLogical = true(size(DGpositions));
DGnFixed = 0;

Dx = deformation_mode;
Dx = zeros(size(deformation_mode));


DGDx = obj.CGxToDGx(Dx);
obj.SetCurrentDGState(DGDx);

v = zeros(length(DGDx),1);

u = [DGDx; v];

K = obj.DGInterfaceStiffnessMatrix;
force = obj.DGInterfaceElasticForce;
b = obj.DGb;

ep = 1;
rng('shuffle')
dx = normc(rand(2*DGN,1));
dx(1) = 0;
dx(3:end) = 0;

obj.SetCurrentDGState(ep*dx + DGDx);
% K_new = obj.DGInterfaceStiffnessMatrix;
% force_new = obj.DGInterfaceElasticForce;
b_new = obj.DGb;
dif = b_new - b;
1;
% disp('max(df-Kdx)');
% disp(max((force_new-force)/ep-K*dx));
% disp(max((force_new-force)/ep-df));
% disp(max(df+K*dx))

end