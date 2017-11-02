function two_elem_eig
%% set parameters, flags, and meshes
close all

rerun_flag = true;
save_state = true;
draw = true;

dt = 1/120;
tsteps = 120*1;

fs = filesep;

mesh_shape = 'two_elem';
simulation_type = 'DGIP';

% set the DG flag base on simulation type
switch simulation_type(1:2)
    case 'DG'
        isDG = true;
    otherwise
        isDG = false;
end

DGeta = 1;
constraints = 1; % types of constraint
% 1: free

Y = 100; % Young's modululs
P = 0.48; % Poisson ratio
rho = 1; % density
a = 0.0; % rayleigh damping
b = 0.00;
material = 'linear'; % choice: 'linear', 'neo-hookean'

axis_box = [-1 1.5 -0.5 1];
elem = [1 2 3; 2 4 3];
nodeM = [0 0; 1 0; 0 1; 1 1]/2;

%%


% % fix the orientation
% elem(:,[1 3]) = elem(:,[3 1]);

N = size(nodeM,1);
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
ne = 6;
% [v,d] = eigs(K,M,ne,'sm');
[v,d] = eig(full(K),full(M));
ind = d > 1e-5;
sort(d(ind));
display('eigenvalue of the stiffness matrix')
display(sort(diag(d)))

if isDG
    % the eigenvalues of the element part
    K_element = obj.DGElementStiffnessMatrix;
    M = obj.DGM;
    [v,d] = eig(full(K_element),full(M));
    %     ind = d > 1e-5;
    %     sort(d(ind))
    display('eigenvalue of the element stiffness matrix')
    display(sort(diag(d)))
    
    
    % the eigenvalues of the glueing part
    K_glue = obj.DGInterfaceStiffnessMatrix;
    M = obj.DGM;
    [v,d] = eig(full(K_glue),full(M));
    %     ind = d > 1e-5;
    %     sort(d(ind))
    display('eigenvalue of the interface stiffness matrix')
    display(sort(diag(d)))
    
end