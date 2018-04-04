function two_elem_eig_script
%% set parameters, flags, and meshes
close all

rerun_flag = true;
save_state = true;
draw = true;

dt = 1/120;
tsteps = 120*1;

fs = filesep;

mesh_shape = 'two_elem';
simulation_type = 'DGBZ';



% set the DG flag base on simulation type
switch simulation_type(1:2)
    case 'DG'
        isDG = true;
        switch simulation_type(3:4)
            case 'BZ'
                isIP = false;
            otherwise
                isIP = true;
        end
    otherwise
        isDG = false;
end
DGeta_list = [1e1 1e0 1e-1 1e-2 1e-3];
DGeta_list = logspace(-2,2);

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
d_list = zeros(length(DGeta_list),1);

for i_DGeta = 1:length(DGeta_list)
    DGeta = DGeta_list(i_DGeta);
    % fix the orientation
%     elem(:,[1 3]) = elem(:,[3 1]);
    
    N = size(nodeM,1);
    % construct triangular mesh object
    if isDG
        obj = elasticDGTriObj(nodeM, elem);
        obj.eta = DGeta;
        switch simulation_type(3:4)
            case 'BZ'
                obj.DGIP = false;
            otherwise
                obj.DGIP = true;
        end
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
    rng(3); % to fix the random seed for eigs
    [v,d] = eigs(K,M,ne,1);
    ind = abs(d) > 1e-5;
    D = sort(d(ind));
    d_list(i_DGeta) = D(1);
%     clear K M v d ne ind D
end

hf = openfig(['result_data' fs 'figure' fs 'two_elem_eig_DG_w_CG_ref']);
hl = findobj('DisplayName','DGIP');
hl.YData = real(d_list');
hl.XData = DGeta_list;
hl.DisplayName = simulation_type;
hl = findobj('DisplayName','CG reference');
hl.XData(1) = min(DGeta_list);
hl.XData(2) = max(DGeta_list);
savefig(hf,['result_data' fs 'figure' fs 'two_elem_eig_DGBZ_w_CG_ref']);

end