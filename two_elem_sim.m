function two_elem_sim
%% set parameters, flags, and meshes
close all

rerun_flag = true;
save_state = true;
draw = true;

dt = 1/120;
tsteps = 120*1;

fs = filesep;

mesh_shape = 'two_elem';
simulation_type = 'CG';

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

DGeta = 1e3;
solver = 'SI';
constraints = 1; % types of constraint
% 1: free
deformation_scale_factor = 10;
deformation_mode_number = 1;

Y = 100; % Young's modululs
P = 0.45; % Poisson ratio
rho = 1; % density
a = 0.0; % rayleigh damping
b = 0.00;
material = 'linear'; % choice: 'linear', 'neo-hookean'

axis_box = [-1 1.5 -0.5 1];

elem = [1 2 3; 2 4 3];
nodeM = [0 0; 1 0; 0 1; 1 1]/2;
%%
% fix the orientation
% elem(:,[1 3]) = elem(:,[3 1]);

N = size(nodeM,1);

dirname = sprintf('sim_data%c%s_%s', fs, material, mesh_shape);
if (exist([dirname fs 'data.mat'], 'file') ~= 2) || rerun_flag
    
    % first construct the CG object for eigen decomps
    obj = elasticTriObj(nodeM, elem);
    switch material
        case 'linear'
            obj.SetMaterial( Y, P, rho, 2, a, b); % set the tri to linear
            %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
        case 'neo-hookean'
            obj.SetMaterial( Y, P, rho, 1, a, b); % set the tri to neo-hookean
    end
    %
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
    
    mkdir(dirname);
    save([dirname fs 'data.mat'], 'obj','V','D'); % storing eigen decomp
else
    load([dirname fs 'data.mat']);
    %     load([filename '.mat'], );
end

% if it is DG, construct triangular mesh object to overwrite the CG
% object
if isDG
    obj = elasticDGTriObj(nodeM, elem);
    obj.eta = DGeta;
    if ~isIP
        obj.DGIP = false;
    end
    switch material
        case 'linear'
            obj.SetMaterial( Y, P, rho, 2, a, b); % set the tri to linear
            %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
        case 'neo-hookean'
            obj.SetMaterial( Y, P, rho, 1, a, b); % set the tri to neo-hookean
    end
end

ha = obj.init_vis;



% deformation for the initial condition

deformation_mode = V(:,3 + deformation_mode_number)/3/deformation_scale_factor;

Dx = deformation_mode;

if isDG
    Dx = obj.CGxToDGx(Dx);
end
obj.SetCurrentState(Dx);
obj.simple_vis(obj.vis_handle);


v = zeros(length(Dx),1);
axis(axis_box)
axis equal
xlim = ha.XLim;
ylim = ha.YLim;
u = [Dx; v];

if save_state && draw
    if isDG
        simdir = strcat(dirname,fs,solver,'_',simulation_type,'_Y',num2str(Y),'_P',num2str(P),'_dt',num2str(dt),'_eta',num2str(DGeta));
    else
        simdir = strcat(dirname,fs,solver,'_',simulation_type,'_Y',num2str(Y),'_P',num2str(P),'_dt',num2str(dt));
    end
    mkdir(simdir);
    vidname = strcat(simdir,fs,'video.avi');
    vid = VideoWriter(vidname);
    vid.FrameRate = 60;
    open(vid);
end

% rate to draw the scene
sim_rate = round(1/dt);
draw_rate = round(sim_rate/vid.FrameRate);

trajectory = zeros(size(u,1),tsteps);
for ti = 1:tsteps
    trajectory(:,ti) = u;
    
    switch solver
        case 'IM'
            u = ImplicitMid(dt, u, obj);
        case 'SI'
            u = SemiImplicit(dt, u, obj);
        case 'SIIMEX'
            u = SemiImplicitIMEX(dt, u, obj);
        case 'ERE'
            u = ERE(dt, u, obj);
    end
    if(draw)
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)
            axis(axis_box)
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


% fname = [filename '_trajectory.mat'];
% save(fname)

save([simdir fs 'trajectory.mat']);



end