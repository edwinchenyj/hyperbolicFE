function rect_simulation

clear all
close all

draw = true;
rerun_flag = true;
save_state = true;
draw = true;

dt = 1/120;
tsteps = 120*3;

fs = filesep;

mesh_shape = 'rect';
maxA = 0.01;
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


DGeta = 1e1;
solver = 'ERE';

Y = 100000; % Young's modululs
P = 0.48; % Poisson ratio
rho = 1000; % density
a = 0.0; % rayleigh damping
b = 0.00;
material = 'neo-hookean'; % choice: 'linear', 'neo-hookean'


axis_box = [-0.5 .5 -3 1];


meshname = sprintf('mesh_data%c%s_maxA_%.d',fs,mesh_shape, maxA);

if exist([meshname '.mat'], 'file') ~= 2
    disp('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end
nodeM = nodeM(:,[2 1]);
elem = elem(:,[1 3 2]);
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
    obj.gravity_on = true;
    obj.calculateGravity;
    %
    Dx = 0*rand(2*N,1); % displacement field. set zero for rest state
    obj.SetCurrentState(Dx);
    mkdir(dirname);
    save([dirname fs 'data.mat'], 'obj'); % storing eigen decomp
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
    obj.gravity_on = true;
    obj.calculateGravity;
end

ha = obj.init_vis;

if ~isDG
    indLogical = true(size(obj.X));
Xind_top = (abs(nodeM(:,2)-max(nodeM(:,2))) < 1e-6);
nFixed = sum(Xind_top);
ind_fix = reshape(transpose(repmat(Xind_top,1,2)),[],1); % logical index for total position vector

indLogical(ind_fix) = false;
end

Dx = zeros(size(obj.X));

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
            u = ImplicitMid(dt, u, obj,~indLogical);
        case 'SI'
            u = SemiImplicit(dt, u, obj,~indLogical);
        case 'SIIMEX'
            u = SemiImplicitIMEX(dt, u, obj, ~indLogical);
        case 'ERE'
            u = ERE(dt, u, obj, ~indLogical);
        case 'BE'
            u = BackwardEuler(dt, u, obj, ~indLogical);
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