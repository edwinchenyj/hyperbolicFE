fname = 'simulation3D_deform_full_SI_not_rescaled_bar3_hete';

parameter_lists = {...
    'save_state',true,...
    'draw',true,...
    'dt',[1/100],...
    'T',[2],...
    'mesh_shape',{'small_bar'},...
    'simulation_type',{'CG'},...
    'solver',{'SI'},...
    'h',[0.48,.3,.2],...
    'Y',2e5,...
    'P',0.45,...
    'rho',1e3,...
    'a',0,...
    'b',0,...
    'material_type',{'neo-hookean'},...
    'gravity',{'on'},...
    'eig_modes',[6,20]};

parse_parameter(fname,[],parameter_lists)