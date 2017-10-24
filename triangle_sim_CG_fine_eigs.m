function triangle_sim_CG_fine_eigs
%% set parameters, flags, and meshes
close all

% draw = false;
rerun_flag = true;
save_state = false;
draw = true;
% test_mode = true;

fs = filesep;

mesh_shape = 'triangle';
maxA = 0.001;
simulation_type = 'CG';

solver = 'IM';
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
ela_material = 'linear'; % choice: 'linear', 'neo-hookean'

axis_box = [-1 1.5 -0.5 1];

meshname = sprintf('mesh_data%c%s_maxA_%.d',fs,mesh_shape, maxA);

if exist([meshname '.mat'], 'file') ~= 2
    disp('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end


%%

N = size(nodeM,1);

dirname = sprintf('sim_data%c%s_%s_maxA_%.d', fs, ela_material, mesh_shape, maxA);
if (exist([dirname fs 'data.mat'], 'file') ~= 2) || rerun_flag
    
    % construct triangular mesh object
        obj = elasticTriObj(nodeM, elem);
    switch ela_material
        case 'linear'
            obj.SetMaterial( Y, P, rho, 1:size(elem,1), 2); % set the tri to linear
            %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
        case 'neo-hookean'
            obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
    end
    %
    Dx = 0*rand(2*N,1); % displacement field. set zero for rest state
    obj.SetCurrentState(Dx);
    %
    M = obj.M;
    K = obj.StiffnessMatrix;
    
    [V,D] = eig(full(K),full(M));
    [low_eig, permutation_indices] = sort(diag(D));
    V = -V(:,permutation_indices);

    
    mkdir(dirname);
    save([dirname fs 'data.mat'], 'obj','V','low_eig'); % storing eigen decomp
else
    load([dirname fs 'data.mat']);
end


end