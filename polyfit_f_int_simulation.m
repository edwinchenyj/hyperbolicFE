function name = polyfit_f_int_simulation(varargin)
close all

draw = true;
rerun_flag = true;
save_state = true;

dt = 1/100;
T = 2;
tsteps = 100*T;

fs = filesep;

mesh_shape = 'circle';
maxA = 0.01;
simulation_type = 'CG';



solver = 'ERE';

Y = 100000; % Young's modululs
P = 0.48; % Poisson ratio
rho = 1000; % density
a = 0.0; % rayleigh damping
b = 0.00;
material_type = 'neo-hookean'; % choice: 'linear', 'neo-hookean'

axis_box = [-0.5 .5 -1.5 1];

deformation_scale = 10; % scale for the initial deformation

gravity = 'off';
mode = 1;

is_polyfit = false;
is_DAC = false;

% parse input
i_arg = 1;
while (i_arg <= nargin)
    switch varargin{i_arg}
        case 'draw'
            i_arg = i_arg + 1;
            draw = varargin{i_arg};
        case 'rerun_flag'
            i_arg = i_arg + 1;
            rerun_flag = varargin{i_arg};
        case 'save_state'
            i_arg = i_arg + 1;
            save_state = varargin{i_arg};
        case 'dt'
            i_arg = i_arg + 1;
            dt = varargin{i_arg};
        case 'T'
            i_arg = i_arg + 1;
            T = varargin{i_arg};
        case 'mesh_shape'
            i_arg = i_arg + 1;
            mesh_shape = varargin{i_arg};
        case 'maxA'
            i_arg = i_arg + 1;
            maxA = varargin{i_arg};
        case 'simulation_type'
            i_arg = i_arg + 1;
            simulation_type = varargin{i_arg};
        case 'DGeta'
            i_arg = i_arg + 1;
            DGeta = varargin{i_arg};
        case 'solver'
            i_arg = i_arg + 1;
            solver = varargin{i_arg};
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
        case 'material'
            i_arg = i_arg + 1;
            material_type = varargin{i_arg};
        case 'axis_box'
            i_arg = i_arg + 1;
            axis_box = varargin{i_arg};
        case 'constraint'
            i_arg = i_arg + 1;
            constraint = varargin{i_arg};
        case 'mode'
            i_arg = i_arg + 1;
            mode = varargin{i_arg};
        case 'gravity'
            i_arg = i_arg + 1;
            gravity = varargin{i_arg};
        case 'deformation_scale'
            i_arg = i_arg + 1;
            deformation_scale = varargin{i_arg};
        case 'polyfit_modes'
            i_arg = i_arg + 1;
            polyfit_modes = varargin{i_arg};
        case 'is_polyfit'
            i_arg = i_arg + 1;
            is_polyfit = varargin{i_arg};
            
    end
    i_arg = i_arg + 1;
end
tsteps = T/dt;
% set the DG flag base on simulation type
switch simulation_type(1:2)
    case 'DG'
        isDG = true;
        %         DGeta = 1e1;
        switch simulation_type(3:4)
            case 'BZ'
                isIP = false;
            otherwise
                isIP = true;
        end
    otherwise
        isDG = false;
end

meshname = sprintf('mesh_data%c%s_maxA_%.d',fs,mesh_shape, maxA);

if exist([meshname '.mat'], 'file') ~= 2
    disp('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end
% nodeM = nodeM(:,[2 1]);
% elem = elem(:,[1 3 2]);
N = size(nodeM,1);


dirname = sprintf('sim_data%c%s_%s', fs, material_type, mesh_shape);

% first construct the CG object for eigen decomps
obj = elasticTriObj(nodeM, elem);
switch material_type
    case 'linear'
        obj.SetMaterial( Y, P, rho, 2, a, b); % set the tri to linear
        %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
    case 'neo-hookean'
        obj.SetMaterial( Y, P, rho, 1, a, b); % set the tri to neo-hookean
end

switch gravity
    case 'on'
        obj.gravity_on = true;
        obj.calculateGravity;
end
%
Dx = 0*rand(2*N,1); % displacement field. set zero for rest state
obj.SetCurrentState(Dx);
mkdir(dirname);

if (exist([dirname fs 'mode_data.mat'], 'file') ~= 2) || rerun_flag
    %
    M = obj.M;
    K = obj.StiffnessMatrix;
    %         K = K(~ind_fix,~ind_fix); % extract the non-fixed part
    
    [V,D] = eig(full(K),full(M));
    [low_eig, permutation_indices] = sort(diag(D));
    V = -V(:,permutation_indices);
    %     firstMode = V(:,4)/2;
    
    save([dirname fs 'mode_data.mat'],'V','D'); % storing eigen decomp
else
    load([dirname fs 'mode_data.mat'],'V','D');
end

% if it is DG, construct triangular mesh object to overwrite the CG
% object
if isDG
    obj = elasticDGTriObj(nodeM, elem);
    obj.eta = DGeta;
    if ~isIP
        obj.DGIP = false;
    end
    switch material_type
        case 'linear'
            obj.SetMaterial( Y, P, rho, 2, a, b); % set the tri to linear
            %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
        case 'neo-hookean'
            obj.SetMaterial( Y, P, rho, 1, a, b); % set the tri to neo-hookean
    end
    %     obj.gravity_on = true;
    %     obj.calculateGravity;
end

ha = obj.init_vis;

if ~isDG
    switch constraint
        case 'right'
            indLogical = true(size(obj.X));
            Xind_side = (abs(nodeM(:,1)-max(nodeM(:,1))) < 1e-6);
            nFixed = sum(Xind_side);
            ind_fix = reshape(transpose(repmat(Xind_side,1,2)),[],1); % logical index for total position vector
            
            indLogical(ind_fix) = false;
            
            Dx = zeros(size(obj.X));
        case 'top'
            indLogical = true(size(obj.X));
            Xind_side = (abs(nodeM(:,2)-max(nodeM(:,2))) < 1e-6);
            nFixed = sum(Xind_side);
            ind_fix = reshape(transpose(repmat(Xind_side,1,2)),[],1); % logical index for total position vector
            
            indLogical(ind_fix) = false;
            
            Dx = zeros(size(obj.X));
        case 'none'
            indLogical = true(size(obj.X));
            
            deformation_mode = V(:,3 + mode)/deformation_scale;
            
            Dx = deformation_mode;
            indLogical = true(size(obj.X));
    end
    
else
    %     indLogical = logical(obj.CGxToDGx(indLogical));
    switch constraint
        case 'none'
            indLogical = true(size(obj.X));
            
            deformation_mode = V(:,3 + mode)/deformation_scale;
            
            Dx = deformation_mode;
            %             indLogical = true(size(obj.X));
    end
    Dx = obj.CGxToDGx(Dx);
end

if (exist([dirname fs 'constrained_mode_data.mat'], 'file') ~= 2) || rerun_flag
    %
    M = obj.M;
    K = obj.StiffnessMatrix;
    M = M(indLogical,indLogical);
    K = K(indLogical,indLogical);
    %         K = K(~ind_fix,~ind_fix); % extract the non-fixed part
    
    [V,D] = eig(full(K),full(M));
    [constrained_low_eig, permutation_indices] = sort(diag(D));
    V = -V(:,permutation_indices);
    %     firstMode = V(:,4)/2;
    
    save([dirname fs 'constrained_mode_data.mat'],'V','D','constrained_low_eig'); % storing eigen decomp
else
    load([dirname fs 'constrained_mode_data.mat'],'V','D','constrained_low_eig');
end

if is_polyfit
    obj.is_polyfit = true;
elseif is_DAC
    obj.is_DAC = true;
end

%% process the high res mesh to extract the polynomial
% set the DG flag base on simulation type
maxA_fine = 0.001;

meshname = sprintf('mesh_data%c%s_maxA_%.d',fs,mesh_shape, maxA_fine);

if exist([meshname '.mat'], 'file') ~= 2
    disp('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end
% nodeM = nodeM(:,[2 1]);
% elem = elem(:,[1 3 2]);
N = size(nodeM,1);


dirname = sprintf('sim_data%c%s_%s', fs, material_type, mesh_shape);

% first construct the CG object for eigen decomps
fine_obj = elasticTriObj(nodeM, elem);
switch material_type
    case 'linear'
        fine_obj.SetMaterial( Y, P, rho, 2, a, b); % set the tri to linear
        %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
    case 'neo-hookean'
        fine_obj.SetMaterial( Y, P, rho, 1, a, b); % set the tri to neo-hookean
end

switch gravity
    case 'on'
        fine_obj.gravity_on = true;
        fine_obj.calculateGravity;
end
%
fine_Dx = 0*rand(2*N,1); % displacement field. set zero for rest state
fine_obj.SetCurrentState(fine_Dx);
mkdir(dirname);

if (exist([dirname fs 'mode_data.mat'], 'file') ~= 2) || rerun_flag
    %
    M = fine_obj.M;
    K = fine_obj.StiffnessMatrix;
    %         K = K(~ind_fix,~ind_fix); % extract the non-fixed part
    
    [V,D] = eig(full(K),full(M));
    [low_eig, permutation_indices] = sort(diag(D));
    V = -V(:,permutation_indices);
    %     firstMode = V(:,4)/2;
    
    save([dirname fs 'mode_data.mat'],'V','D'); % storing eigen decomp
else
    load([dirname fs 'mode_data.mat'],'V','D');
end


ha = fine_obj.init_vis;

if ~isDG
    switch constraint
        case 'right'
            fine_indLogical = true(size(fine_obj.X));
            Xind_side = (abs(nodeM(:,1)-max(nodeM(:,1))) < 1e-6);
            nFixed = sum(Xind_side);
            ind_fix = reshape(transpose(repmat(Xind_side,1,2)),[],1); % logical index for total position vector
            
            fine_indLogical(ind_fix) = false;
            
            fine_Dx = zeros(size(fine_obj.X));
        case 'top'
            fine_indLogical = true(size(fine_obj.X));
            Xind_side = (abs(nodeM(:,2)-max(nodeM(:,2))) < 1e-6);
            nFixed = sum(Xind_side);
            ind_fix = reshape(transpose(repmat(Xind_side,1,2)),[],1); % logical index for total position vector
            
            fine_indLogical(ind_fix) = false;
            
            fine_Dx = zeros(size(fine_obj.X));
        case 'none'
            fine_indLogical = true(size(fine_obj.X));
            
            deformation_mode = V(:,3 + mode)/deformation_scale;
            
            fine_Dx = deformation_mode;
            fine_indLogical = true(size(fine_obj.X));
    end
    
else
    %     indLogical = logical(obj.CGxToDGx(indLogical));
    switch constraint
        case 'none'
            fine_indLogical = true(size(fine_obj.X));
            
            deformation_mode = V(:,3 + mode)/deformation_scale;
            
            fine_Dx = deformation_mode;
            %             indLogical = true(size(obj.X));
    end
    fine_Dx = fine_obj.CGxToDGx(fine_Dx);
end


if (exist([dirname fs 'constrained_mode_data.mat'], 'file') ~= 2) || rerun_flag
    %
    M = fine_obj.M;
    K = fine_obj.StiffnessMatrix;
    M = M(fine_indLogical,fine_indLogical);
    K = K(fine_indLogical,fine_indLogical);
    %         K = K(~ind_fix,~ind_fix); % extract the non-fixed part
    
    [V,D] = eig(full(K),full(M));
    [constrained_low_eig_fine, permutation_indices] = sort(diag(D));
    V = -V(:,permutation_indices);
    %     firstMode = V(:,4)/2;
    
    save([dirname fs 'constrained_mode_data.mat'],'V','D','constrained_low_eig_fine'); % storing eigen decomp
else
    load([dirname fs 'constrained_mode_data.mat'],'V','D','constrained_low_eig_fine');
end
obj.polyfit_p = polyfit([0;constrained_low_eig(1:polyfit_modes)],[0;constrained_low_eig_fine(1:polyfit_modes)],polyfit_modes);
obj.polyfit_p(end) = 0;
% obj.polyfit_p = polyfit([constrained_low_eig(1:polyfit_modes)],[constrained_low_eig_fine(1:polyfit_modes)],polyfit_modes-1);
%%
obj.SetCurrentState(Dx);
% obj.simple_vis(obj.vis_handle);


v = zeros(length(Dx),1);
axis(ha,axis_box)
axis equal
xlim = ha.XLim;
ylim = ha.YLim;
u = [Dx; v];
if sum(Dx) == 0
    obj.from_rest = true;
    obj.polyfit_force_approx = zeros(sum(indLogical),1);
end

if save_state && draw
    if isDG
        simdir = strcat(dirname,fs,solver,'_',...
            'polyfit_f_int_',...
            simulation_type,...
            '_constraint_',constraint,...
            '_maxA',num2str(maxA),...
            '_Y',num2str(Y),...
            '_P',num2str(P),...
            '_rho',num2str(rho),...
            '_a',num2str(a),...
            '_b',num2str(b),...
            '_dt',num2str(dt),...
            '_eta',num2str(DGeta),...
            '_def-scl',num2str(deformation_scale),...
            'polyfit_modes',num2str(polyfit_modes));
    else
        simdir = strcat(dirname,fs,solver,'_',simulation_type,...
            'polyfit_f_int_',...
            '_constraint_',constraint,...
            '_maxA',num2str(maxA),...
            '_Y',num2str(Y),...
            '_P',num2str(P),...
            '_rho',num2str(rho),...
            '_a',num2str(a),...
            '_b',num2str(b),...
            '_dt',num2str(dt),...
            '_def-scl',num2str(deformation_scale),...
            'polyfit_modes',num2str(polyfit_modes));
        
    end
    mkdir(simdir);
    vidname = strcat(simdir,fs,'video.avi');
    vid = VideoWriter(vidname);
    vid.FrameRate = 50;
    open(vid);
end

elastic_energy = zeros(1,tsteps);
kinetic_energy = zeros(1,tsteps);
gravitational_potential = zeros(1,tsteps);

% rate to draw the scene
sim_rate = round(1/dt);
draw_rate = round(sim_rate/vid.FrameRate);

if (exist([simdir fs 'trajectory.mat'], 'file') ~= 2)
    trajectory = zeros(size(u,1),tsteps);
    for ti = 1:tsteps
        trajectory(:,ti) = u;
        elastic_energy(ti) = obj.ElasticEnergy;
        if ~is_polyfit
            switch solver
                case 'IM'
                    u = ImplicitMid(dt, u, obj,~indLogical);
                case 'SI'
                    u = SemiImplicit(dt, u, obj,~indLogical);
                case 'SIModified'
                    u = SemiImplicitModified(dt, u, obj,~indLogical);
                case 'SIIMEX'
                    u = SemiImplicitIMEX(dt, u, obj, ~indLogical);
                case 'SIIMEXModified'
                    u = SemiImplicitIMEXModified(dt, u, obj, ~indLogical);
                case 'SIRK2IMEX'
                    u = SemiImplicitRK2IMEX(dt, u, obj, ~indLogical);
                case 'SIRK4IMEX'
                    u = SemiImplicitRK4IMEX(dt, u, obj, ~indLogical);
                case 'SIEXPINTIMEX'
                    u = SemiImplicitEXPINTIMEX(dt, u, obj, ~indLogical);
                case 'ERE'
                    u = ERE(dt, u, obj, ~indLogical);
                case 'BE'
                    u = BackwardEuler(dt, u, obj, ~indLogical);
                case 'IMIMEX'
                    u = IMIMEX(dt, u, obj, ~indLogical);
                case 'BEIMEX'
                    u = BEIMEX(dt, u, obj, ~indLogical);
            end
            
        else
            switch solver
                case 'IM'
                    u = polyfit_f_int_ImplicitMid(dt, u, obj,~indLogical);
                case 'SI'
                    u = polyfit_f_int_SemiImplicit(dt, u, obj,~indLogical);
                case 'ERE'
                    u = polyfit_f_int_ERE(dt, u, obj, ~indLogical);
                case 'BE'
                    u = polyfit_f_int_BackwardEuler(dt, u, obj, ~indLogical);
            end
                            
        end
        if(draw)
            if or(mod(ti, draw_rate) == 1, draw_rate == 1)
                axis(ha,axis_box)
                cla
                obj.simple_vis(obj.vis_handle);
                obj.vis_handle.XLim = xlim;
                obj.vis_handle.YLim = ylim;
                if save_state
                    frame = getframe(gca);
                    writeVideo(vid,frame);
                end
            end
            
        end
    end
    
    % plot the energy
    for ti = 1:tsteps
        u = trajectory(:,ti);
        x = u(1:end/2);
        v = u(end/2+1:end);
        gravitational_potential(ti) = -(x+obj.X)' * obj.M * obj.externalGravity;
        kinetic_energy(ti) = 1/2 * v' * obj.M * v;
    end
    cla reset;
    plot(1:tsteps,elastic_energy,1:tsteps,kinetic_energy,1:tsteps,gravitational_potential,...
        1:tsteps,elastic_energy+kinetic_energy+gravitational_potential)
    print([simdir fs 'energy'],'-dpng')
    save([simdir fs 'trajectory.mat']);
else
    disp([material_type simdir fs 'trajectory.mat'])
    disp('already exist. delete it if you want to re-run')
end

name = [simdir fs 'trajectory.mat'];
close(vid);
end