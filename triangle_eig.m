function triangle_eig
%% set parameters, flags, and meshes
close all

rerun_flag = true;
save_state = true;
draw = true;

dt = 1/120;
tsteps = 120*1;

fs = filesep;

mesh_shape = 'triangle';
maxA = 0.1;
simulation_type = 'DGBZ';

% set the DG flag base on simulation type
switch simulation_type(1:2)
    case 'DG'
        isDG = true;
    otherwise
        isDG = false;
end

DGeta = 1e-1;
solver = 'SI';
constraints = 1; % types of constraint
% 1: free

deformation_mode_number = 1;
switch maxA
    case 0.1
        deformation_scale_factor = 2;
    case 0.01
        deformation_scale_factor = -2; % there is a sign change when maxA = 0.01, 0.001
    case 0.001
        deformation_scale_factor = -2; % there is a sign change when maxA = 0.01, 0.001
end
Y = 100; % Young's modululs
P = 0.48; % Poisson ratio
rho = 1; % density
a = 0.0; % rayleigh damping
b = 0.00;
material = 'linear'; % choice: 'linear', 'neo-hookean'

axis_box = [-1 1.5 -0.5 1];

meshname = sprintf('mesh_data%c%s_maxA_%.d',fs,mesh_shape, maxA);

if exist([meshname '.mat'], 'file') ~= 2
    disp('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
end

% options for fsolve in implicit solver
options = optimoptions('fsolve','TolFun',1.e-9,'TolX',1.e-9,'Display','final');
% options = optimoptions('fsolve','Algorithm','levenberg-marquardt');

%%


% fix the orientation
% [nodeM, elem] = fixmesh(nodeM, elem);
elem(:,[1 3]) = elem(:,[3 1]);

N = size(nodeM,1);

dirname = sprintf('sim_data%c%s_%s_maxA_%.d', fs, material, mesh_shape, maxA);

% construct triangular mesh object
if isDG
    obj = elasticDGTriObj(nodeM, elem);
    obj.eta = DGeta;
else
    obj = elasticTriObj(nodeM, elem);
end

switch material
    case 'linear'
        obj.SetMaterial( Y, P, rho, 1:size(elem,1), 2); % set the tri to linear
        %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
    case 'neo-hookean'
        obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
end


% deformation for the initial condition

Dx = zeros(length(obj.X),1);

if isDG
    N = size(obj.DGnodeM,1);
    node = obj.DGnode;
else
    
    node = obj.node;
end

positionsM = node';
positionsM = positionsM(:);
positions = positionsM;
externalGravity  = zeros(size(positions));

if isDG
    indLogical = true(size(positions)); % TODO: need to map the indLogical to the indLogical for DG
    nFixed = 0;
    
    Dx = obj.CGxToDGx(Dx);
    obj.SetCurrentDGState(Dx);
    
    K = obj.DGStiffnessMatrix;
    M = obj.DGM;
else
    indLogical = true(size(positions)); % TODO: need to map the indLogical to the indLogical for DG
    nFixed = 0;
    obj.SetCurrentState(Dx);
    
    K = obj.StiffnessMatrix;
    M = obj.M;
    
end
ne = 10;
[v,d] = eigs(K,M,ne,'sm');
sort(diag(d))
end