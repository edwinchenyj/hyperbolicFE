function name = simulation3D_motions_bar(varargin)
close all
clf

draw = true;
rerun_flag = true;
save_state = true;

dt = 1/100;
T = 4;
tsteps = 100*T;

fs = filesep;

mesh_shape = 'small_bar';
h = 0.2;
simulation_type = 'CG';



solver = 'SI';

Y = 100000; % Young's modululs
P = 0.48; % Poisson ratio
rho = 1000; % density
a = 0.0; % rayleigh damping
b = 0.00;
material_type = 'neo-hookean'; % choice: 'linear', 'neo-hookean'

axis_box = [-0.5 .5 -1.5 1];

deformation_scale = 10; % scale for the initial deformation

gravity = 'off';
mode = 'n';
constraint = 'n';
deformation_scale = 'n';
DGeta = 'n';

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
        case 'h'
            i_arg = i_arg + 1;
            h = varargin{i_arg};
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
            
    end
    i_arg = i_arg + 1;
end

tsteps = T/dt;


meshname = sprintf('mesh_data%c%s_h_%.d',fs,mesh_shape, h);

if exist([meshname '.mat'], 'file') ~= 2
    disp('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end
nodeM = nodeM(:,[3 2 1]);
elem = elem(:,[1 2 3 4]);
N = size(nodeM,1);

%% load Dx from free vibration bar
% dirname = sprintf('sim_data%c%s_%s_%s', fs, 'simulation3D_free_bar',material_type, mesh_shape);
% Y_vibration = 1e4;
% P_vibration = 0.45;
% rho_vibration = 1000;
% a_vibration = 0.0;
% b_vibration = 0.00;
% dt_vibration = 1/100;
% simdir = strcat(dirname,fs,'BE','_',simulation_type,...
%     '_constraint_',constraint,...
%     '_h',num2str(h),...
%     '_Y',num2str(Y_vibration),...
%     '_P',num2str(P_vibration),...
%     '_rho',num2str(rho_vibration),...
%     '_a',num2str(a_vibration),...
%     '_b',num2str(b_vibration),...
%     '_dt',num2str(dt_vibration),...
%     '_def-scl',num2str(deformation_scale));
% 
% 
% 
% vibration = load([simdir fs 'deformation.mat'], 'Dx'); % storing eigen decomp
% vibration_def = reshape(vibration.Dx,3,[]);

% %% load Dx from dropping bar
% traj_name = ['sim_data\simulation3D_dropping_bar_',material_type,'_small_bar\ERE_CG_constraint_n_h',num2str(h),'_Y100000_P0.45_rho1000_a0.01_b0.01_dt0.01_def-scln\trajectory.mat'];
% 
% dropping_traj = load(traj_name);
% close(gcf)
% dropping_dx = reshape(dropping_traj.trajectory(1:end/2,200),3,[]);
% dropping_dx = dropping_dx([3 2 1],:);
% dropping_dx = dropping_dx(:);
% %% load Dx from releasing bar
% dirname = sprintf('sim_data%c%s_%s_%s', fs, 'simulation3D_releasing_bar',material_type, mesh_shape);
% Y_releasing = 1e3;
% P_releasing = 0.45;
% rho_releasing = 1000;
% a_releasing = 0.1;
% b_releasing = 0.01;
% dt_releasing = 1/100;
% simdir = strcat(dirname,fs,'ERE','_',simulation_type,...
%     '_constraint_',constraint,...
%     '_h',num2str(h),...
%     '_Y',num2str(Y_releasing),...
%     '_P',num2str(P_releasing),...
%     '_rho',num2str(rho_releasing),...
%     '_a',num2str(a_releasing),...
%     '_b',num2str(b_releasing),...
%     '_dt',num2str(dt_releasing),...
%     '_def-scl',num2str(deformation_scale));
% 
% 
% 
% twisting = load([simdir fs 'deformation.mat'], 'Dx','ind','R'); % storing eigen decomp
% 
% % turn twisting def 90 degree
% twisting_dx = reshape(twisting.Dx,3,[]);
% twisting_dx = twisting_dx([3,2,1],:);
% twisting.Dx = twisting_dx(:);
% 
% % put the vibration_def to the sorted order 
% % vibration_def = vibration_def(:,twisting.ind);
% % vibration.Dx = vibration_def(:);
% 
% elem = twisting.R(elem(:,:));

dirname = sprintf('sim_data%c%s_%s_%s', fs, mfilename,material_type, mesh_shape);

mkdir(dirname);
simdir = strcat(dirname,fs,solver,'_',simulation_type,...
    '_constraint_',constraint,...
    '_h',num2str(h),...
    '_Y',num2str(Y),...
    '_P',num2str(P),...
    '_rho',num2str(rho),...
    '_a',num2str(a),...
    '_b',num2str(b),...
    '_dt',num2str(dt),...
    '_def-scl',num2str(deformation_scale));

mkdir(simdir);

%% set fixed points and
% fix nFixed number of highest point

% nodeM = nodeM(twisting.ind,:);
% elem = twisting.R(elem(:,:));

bottom_points = find(abs(nodeM(:,3)-0) < 0.01);
top_points = find(abs(nodeM(:,3)-1/5) < 0.01/5);

%% get logical indices 
indBottom = bottom_points;
indTop = top_points;

nFixed = length(indTop);

N = size(nodeM,1); % number of nodes
indAll = 1:N;
indRemove = indAll([indTop(:)]);
indLogical = logical(ones(3,N));
indLogical(:,indRemove) = logical(0);
indLogical = indLogical(:);

indLeftLogical = logical(zeros(3,N));
indLeftLogical(:,indBottom) = logical(1);
indLeftLogical = indLeftLogical(:);

indRightLogical = logical(zeros(3,N));
indRightLogical(:,indTop) = logical(1);
indRightLogical(:);


obj = elasticTetObj(nodeM, elem);
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
    case 'off'
        obj.gravity_on = false;
        obj.calculateGravity;
end







%
% Dx = 0*rand(obj.Dim*N,1); % displacement field. set zero for rest state
% Dx = loaded.obj.x - loaded.obj.X;
% Dx = twisting.Dx + dropping_dx * 2;
% Dx = twisting.Dx;
% Dx = vibration.Dx;
% temp_Dx = Dx(~indLogical);
% temp_Dx(3:3:end) = max(temp_Dx);
% Dx(~indLogical) = temp_Dx;
obj.SetCurrentState(Dx);

ha = obj.init_vis;
obj.simple_vis(obj.vis_handle);

v = zeros(length(Dx),1);
axis(ha,axis_box)
axis equal
% xlim = ha.XLim;
% ylim = ha.YLim;
xlim = axis_box(1:2);
ylim = axis_box(3:4);
zlim = axis_box(5:6);
u = [Dx; v];
obj.vis_handle.XLim = xlim;
obj.vis_handle.YLim = ylim;
obj.vis_handle.ZLim = zlim;
camP = campos(gca);
camV = camva(gca);

fig = gcf;

fig.Children.Visible = 'off';
fig.Children.Clipping = 'off';
fig.Children.Projection = 'perspective';
fig.Position(1:2) = [1300 50];
fig.Position(3:4) = fig.Position(3:4) *2;
skymap = colormap(gray/4+0.3); %We'll use the colormap "winter" which is defined by matlab.

% Here's a gradient I made up that's a bit more horizon-like. You can
% experiment.
% flipud([logspace(-0.5,0,1000)',logspace(-0.5,0,1000)',logspace(-0.001,0,1000)']);

[skyX,skyY,skyZ] = sphere(200); %create the surface data for a sphere.
% sky = surf(500000*skyX,500000*skyY,500000*skyZ,'LineStyle','none','FaceColor','interp'); %Make a sphere object.

if save_state && draw
    

    vidname = strcat(simdir,fs,'video.avi');
%     vid = VideoWriter(vidname,'Uncompressed AVI');
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

patch_color = [0.8,0.9,0.5];
if (exist([simdir fs 'trajectory.mat'], 'file') ~= 2)
    trajectory = zeros(size(u,1),tsteps);
    for ti = 1:tsteps
        trajectory(:,ti) = u;
        %         elastic_energy(ti) = obj.ElasticEnergy;
        
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
            case 'RK4'
                u = RK4(dt, u, obj, ~indLogical);
        end
        

        if(draw)
            if or(mod(ti, draw_rate) == 1, draw_rate == 1)
                axis(ha,axis_box)
                cla
                
                
                obj.simple_vis(obj.vis_handle,patch_color);
                obj.vis_handle.XLim = xlim;
                obj.vis_handle.YLim = ylim;
                obj.vis_handle.ZLim = zlim;
                
                ah = gca;
                ph = ah.Children;
                ph.FaceLighting = 'flat';
                %                 ph.EdgeLighting = 'gouraud';
                %                 ph.BackFaceLighting = 'lit';
                lighting flat
                material metal
                %                 ph.EdgeColor = [0.9 0.9 1];
                lh = camlight('headlight');
                lh.Style = 'local';
                %                 lh2 = camlight('headlight');
                %                 lh.Position(2) = 3;
                %                 lh.Position(3) = 3;
                map = [0.8 0.9 1];
                %                 colormap(autumn);
                
                shading interp
                %                 ph.EdgeColor = [0.1 0.1 0.1];
                ph.EdgeLighting = 'none';
                ph.EdgeAlpha = 0;
                %                 ph.LineWidth = 0.25;
                %                 ph.AlignVertexCenters = 'on';
                whitebg('white')
                
                %                 origPos = [0,0,0]';
                % %                 campos(origPos);
                %                 origTarget = [0,0,-10000];
                %                 camtarget(origTarget);
                %                 camva(40)
                
                fig.Children.Visible = 'off';
                fig.Children.Clipping = 'off';
                fig.Children.Projection = 'perspective';
                campos([-1.1513   -1.6081    1.4854]);
                camtarget([0.0810   -0.0021    0.0680])
                camva(6.9295);
                scal = 2000;
                offc = 0.1;
                hold on
                sky = surf(scal*(skyX-offc),scal*(skyY-offc),scal*(skyZ-offc),'LineStyle','none','FaceColor','flat','FaceLighting','gouraud'); %Make a sphere object.
                lh.Style = 'local';
                patch_CData = obj.patch_handle.CData;
                colormap(winter);
                % obj.patch_handle.FaceVertexCData = ones(63,3) *0.5;
                obj.patch_handle.FaceVertexCData = repmat([0.9100 0.4100 0.1700],N,1);
                % obj.patch_handle.FaceVertexCData = patch_CData';
                if save_state
                    frame = getframe(fig);
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