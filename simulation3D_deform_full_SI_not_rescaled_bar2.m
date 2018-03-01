function name = simulation3D_deform_full_SI_not_rescaled_bar2(varargin)
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

eig_modes = 6;

solver = 'ERE';

Y = 100000; % Young's modululs
P = 0.48; % Poisson ratio
rho = 1000; % density
a = 0.0; % rayleigh damping
b = 0.00;
material_type = 'neo-hookean'; % choice: 'linear', 'neo-hookean'

axis_box = [-0.5 .5 -1.5 1 -1 1];

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
        case 'reduced_dimension'
            i_arg = i_arg + 1;
            reduced_dimension = varargin{i_arg};
            
    end
    i_arg = i_arg + 1;
end

tsteps = T/dt;

draw = false;
%% process the high res mesh first
h_original = h;
% h = 0.1;
% meshname = sprintf('mesh_data%c%s_h_%.d',fs,mesh_shape, h);
% 
% if exist([meshname '.mat'], 'file') ~= 2
%     disp('mesh does not exist')
%     
% else
%     load([meshname '.mat'], 'nodeM', 'elem');
%     
% end
% % nodeM = nodeM(:,[2 1]);
% % elem = elem(:,[1 3 2]);
% N = size(nodeM,1);
% % 
% % loaddir = sprintf('sim_data%csimulation3D_rotate_bar_%s_%s', fs,material_type, mesh_shape);
% % simdir = strcat(loaddir,fs,'SI','_','CG',...
% %     '_constraint_',constraint,...
% %     '_h',num2str(h),...
% %     '_Y',num2str(1e5),...
% %     '_P',num2str(0.45),...
% %     '_rho',num2str(1000),...
% %     '_a',num2str(0.01),...
% %     '_b',num2str(0.00),...
% %     '_dt',num2str(1/100),...
% %     '_def-scl',num2str(deformation_scale));
% % loaded =load([simdir fs 'trajectory.mat'],'obj');
% 
% % dirname = sprintf('sim_data%c%s_%s_%s', fs, mfilename,material_type, mesh_shape);
% 
% % fix nFixed number of highest point
% % fix nFixed number of highest point
% 
% [~, ind] = sortrows(nodeM,1);
% 
% % put the nodes in order (in z)
% nodeM = nodeM(ind,:);
% 
% % create the rank for the original z
% [~,R] = sort(ind);
% % the new element list
% elem = R(elem(:,:));
% 
% % nodeM = loaded.obj.nodeM;
% % elem = loaded.obj.elem;
% 
% left_points = find(abs(nodeM(:,1)-0) < 0.01);
% right_points = find(abs(nodeM(:,1)-1/5) < 0.01/5);
% 
% % get logical indices for the left and right nodes
% 
% indLeft = left_points;
% indRight = right_points;
% 
% nFixed = length(indRight);
% 
% N = size(nodeM,1); % number of nodes
% indAll = 1:N;
% indRemove = indAll([indRight(:)]);
% indLogical = logical(ones(3,N));
% indLogical(:,indRemove) = logical(0);
% indLogical = indLogical(:);
% 
% indLeftLogical = logical(zeros(3,N));
% indLeftLogical(:,indLeft) = logical(1);
% indLeftLogical = indLeftLogical(:);
% 
% indRightLogical = logical(zeros(3,N));
% indRightLogical(:,indRight) = logical(1);
% indRightLogical(:);
% 
% 
% obj = elasticTetObj(nodeM, elem);
% switch material_type
%     case 'linear'
%         obj.SetMaterial( Y, P, rho, 2, a, b); % set the tri to linear
%         %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
%     case 'neo-hookean'
%         obj.SetMaterial( Y, P, rho, 1, a, b); % set the tri to neo-hookean
%     case 'stvk'
%         obj.SetMaterial( Y, P, rho, 3, a, b); % set the tri to stvk
% end
% 
% 
% Dx = 0*rand(obj.Dim*N,1); % displacement field. set zero for rest state
% obj.SetCurrentState(Dx);
% %
% M = obj.M;
% K = obj.StiffnessMatrix;
% M = M(indLogical,indLogical);
% K = K(indLogical,indLogical);
% %         K = K(~ind_fix,~ind_fix); % extract the non-fixed part
% 
% [V,D] = eigs(K,M,eig_modes,'sm');
% [low_eig, permutation_indices] = sort(diag(D));
% 
% high_res_low_eig = low_eig;

%%
h = h_original;
meshname = sprintf('mesh_data%c%s_h_%.2d',fs,mesh_shape, h);

if exist([meshname '.mat'], 'file') ~= 2
    err('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end
% nodeM = nodeM(:,[3 2 1]);
% elem = elem(:,[1 2 3 4]);

% %% load Dx from dropping bar
% traj_name = ['sim_data' fs 'simulation3D_dropping_bar_',material_type,'_small_bar' fs 'ERE_CG_constraint_n_h',num2str(h),'_Y100000_P0.45_rho1000_a0.01_b0.01_dt0.01_def-scln' fs 'trajectory.mat'];
% 
% dropping_traj = load(traj_name);
% close(gcf)
% dropping_dx = reshape(dropping_traj.trajectory(1:end/2,200),3,[]);
% dropping_dx = dropping_dx([3 2 1],:);
% dropping_dx = dropping_dx(:);

%% load Dx from rotating bar

traj_name = ['sim_data' fs 'simulation3D_rotate_full_SI_rescaled_bar_',...
    'stvk','_small_bar' fs 'SI_CG_constraint_n_h',num2str(h),'_Y10000_P0.45_rho1000_a0_b0_dt0.01_def-scln' fs 'trajectory.mat'];
% 
twisting_traj = load(traj_name);
close(gcf)
twisting_dx = reshape(twisting_traj.trajectory(1:end/2,end),3,[]);
% twisting_dx = twisting_dx([1 2 1],:);
twisting_dx = twisting_dx(:);

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
% twisting = load([simdir fs 'deformation.mat'], 'Dx','ind','R'); % storing eigen decomp
% 
% % turn twisting def 90 degree
% twisting_dx = reshape(twisting.Dx,3,[]);
% twisting_dx = twisting_dx([3,2,1],:);
% twisting.Dx = twisting_dx(:);

% put the vibration_def to the sorted order 
% vibration_def = vibration_def(:,twisting.ind);
% vibration.Dx = vibration_def(:);

% elem = twisting.R(elem(:,:));

%%
dirname = sprintf('sim_data%c%s_%s_%s', fs, [mfilename num2str(eig_modes)],material_type, mesh_shape);

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
    '_def-scl',num2str(deformation_scale),...
    '_redu',num2str(reduced_dimension));

mkdir(simdir);


% nodeM = nodeM(:,[3 2 1]);
elem = elem(:,[1 2 3 4]);
N = size(nodeM,1);

%% set fixed points and
% fix nFixed number of highest point

[~, ind] = sortrows(nodeM,1);

% put the nodes in order (in z)
nodeM = nodeM(ind,:);

% create the rank for the original z
[~,R] = sort(ind);
% the new element list
elem = R(elem(:,:));

bottom_points = find(abs(nodeM(:,1)-0) < 0.01);
top_points = find(abs(nodeM(:,1)-1/5) < 0.01/5);

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

% nodeM = nodeM(:,[3 2 1]);
% elem = elem(:,[1 2 3 4]);

obj = elasticTetObj(nodeM, elem);
switch material_type
    case 'linear'
        obj.SetMaterial( Y, P, rho, 2, a, b); % set the tri to linear
        %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
    case 'neo-hookean'
        obj.SetMaterial( Y, P, rho, 1, a, b); % set the tri to neo-hookean
    case 'stvk'
        obj.SetMaterial( Y, P, rho, 3, a, b); % set the tri to stvk
end

switch gravity
    case 'on'
        obj.gravity_on = true;
        obj.calculateGravity;
    case 'off'
        obj.gravity_on = false;
        obj.calculateGravity;
end


Dx = 0*rand(obj.Dim*N,1); % displacement field. set zero for rest state
obj.SetCurrentState(Dx);

obj.indLogical = indLogical;

% obj.high_res_eig_to_match = high_res_low_eig;
% obj.reduced_dimension = reduced_dimension;
% 
% M = obj.M;
% K = obj.StiffnessMatrix;
% M = M(indLogical,indLogical);
% K = K(indLogical,indLogical);
% 

% [V,D] = eigs(K,M,eig_modes,'sm');
% V = V(:,permutation_indices);
% % if h ~= 0.5
% %     V = -V;
% % end
% [low_eig, permutation_indices] = sort(diag(D));
% obj.eig_ratios = high_res_low_eig./low_eig;
% obj.eig_targets = high_res_low_eig;
% obj.eig_modes = eig_modes;
% Dx_mode = Dx;
% Dx_mode(indLogical) = sum([V(:,1)/20],2);
Dx = twisting_dx;
% Dx = 3.5*dropping_dx;
% Dx = zeros(size(dropping_dx * 2));

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
fig.Position(1:2) = [100 50];
fig.Position(3:4) = fig.Position(3:4) *2;
skymap = colormap(gray/4+0.3); %We'll use the colormap "winter" which is defined by matlab.

[skyX,skyY,skyZ] = sphere(200); %create the surface data for a sphere.
% sky = surf(500000*skyX,500000*skyY,500000*skyZ,'LineStyle','none','FaceColor','interp'); %Make a sphere object.

if save_state && draw
    
    vidname = strcat(simdir,fs,'video.avi');
%     vid = VideoWriter(vidname,'Uncompressed AVI');
    vid = VideoWriter(vidname);
    vid.FrameRate = 50;
    open(vid);
% rate to draw the scene
sim_rate = round(1/dt);
draw_rate = round(sim_rate/vid.FrameRate);

end

elastic_energy = zeros(1,tsteps);
kinetic_energy = zeros(1,tsteps);
gravitational_potential = zeros(1,tsteps);


patch_color = [0.8,0.9,0.5];
if (exist([simdir fs 'trajectory.mat'], 'file') ~= 2)
    trajectory = zeros(size(u,1),tsteps);
    for ti = 1:tsteps
        trajectory(:,ti) = u;
        %         elastic_energy(ti) = obj.ElasticEnergy;
        
        switch solver
            case 'SI'
                u = SemiImplicit(dt, u, obj,~indLogical);
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
if draw
close(vid);
end
end