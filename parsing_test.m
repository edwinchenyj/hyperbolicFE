parameter_lists = {'dt',[1/100],...
    'T',[2],...
    'mesh_shape',{'small_bar'},...
    'simulation_type',{'CG'},...
    'solver',{'SI'},...
    'h',[0.48,0.3,0.2],...
    'Y',2e5,...
    'P',0.45,...
    'rho',1e3,...
    'a',0,...
    'b',0,...
    'material_type',{'neo-hookean'},...
    'gravity',{'on'},...
    'eig_modes',[6,20]};

parse_parameter('func',[],parameter_lists)