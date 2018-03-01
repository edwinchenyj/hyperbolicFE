
function sim_annotation = simulation3D_scene_annotation(varargin)
close all
clf

% parse input

fs = filesep;

i_arg = 1;
while (i_arg <= nargin)
    switch varargin{i_arg}
        
        case 'dt'
            i_arg = i_arg + 1;
            dt = varargin{i_arg};
        case 'T'
            i_arg = i_arg + 1;
            T = varargin{i_arg};
        case 'meshfile'
            i_arg = i_arg + 1;
            meshfile = varargin{i_arg};
        case 'h'
            i_arg = i_arg + 1;
            h = varargin{i_arg};
            %         case 'simulation_type'
            %             i_arg = i_arg + 1;
            %             simulation_type = varargin{i_arg};
        case 'DGeta'
            i_arg = i_arg + 1;
            DGeta = varargin{i_arg};
        case 'solver'
            i_arg = i_arg + 1;
            solver = eval(['@' varargin{i_arg}]);
        case 'Y'
            i_arg = i_arg + 1;
            Y = varargin{i_arg};
        case 'P'
            i_arg = i_arg + 1;
            P = varargin{i_arg};
        case 'rho'
            i_arg = i_arg + 1;
            rho = varargin{i_arg};
        case 'a'
            i_arg = i_arg + 1;
            a = varargin{i_arg};
        case 'b'
            i_arg = i_arg + 1;
            b = varargin{i_arg};
        case 'material_type'
            i_arg = i_arg + 1;
            material_type = varargin{i_arg};
        case 'axis_box'
            i_arg = i_arg + 1;
            axis_box = varargin{i_arg};
        case 'constraint'
            i_arg = i_arg + 1;
            constraint = varargin{i_arg};
        case 'eig_modes'
            i_arg = i_arg + 1;
            eig_modes = varargin{i_arg};
        case 'gravity'
            i_arg = i_arg + 1;
            gravity = varargin{i_arg};
        case 'deformation_scale'
            i_arg = i_arg + 1;
            deformation_scale = varargin{i_arg};
        case 'reduced_dimension'
            i_arg = i_arg + 1;
            reduced_dimension = varargin{i_arg};
        case 'NMSC'
            i_arg = i_arg + 1;
            NMSC = varargin{i_arg};
        case 'fine_meshfile'
            i_arg = i_arg + 1;
            fine_meshfile = varargin{i_arg};
        case 'scene_name'
            i_arg = i_arg + 1;
            scene_name = varargin{i_arg};
            
    end
    i_arg = i_arg + 1;
end

tsteps = T/dt;
%%
scriptdir = script_directory(mfilename);
mkdir(scriptdir);
sim_directory_name = sim_directory(varargin);
simdir = [scriptdir,fs,sim_directory_name];
mkdir(simdir);

str_ind = strfind(sim_directory_name,'__');
for i_str_ind = 1:2:length(str_ind)
    sim_directory_name(str_ind(i_str_ind):str_ind(i_str_ind)+1) =  '=';
end
for i_str_ind = 2:2:length(str_ind)
    sim_directory_name(str_ind(i_str_ind):str_ind(i_str_ind)+1) =  newline;
end
sim_directory_name = strrep(sim_directory_name,'==','=');
sim_directory_name = strrep(sim_directory_name,[newline newline],newline);
sim_annotation = sim_directory_name;

end