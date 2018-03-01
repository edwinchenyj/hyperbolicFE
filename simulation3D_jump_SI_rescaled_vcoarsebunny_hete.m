function name = simulation3D_jump_SI_rescaled_vcoarsebunny(varargin)
close all
clf

draw = true;
rerun_flag = true;
save_state = true;

dt = 1/100;
T = 4;
tsteps = 100*T;

fs = filesep;

eig_modes = 20;

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
            
    end
    i_arg = i_arg + 1;
end

tsteps = T/dt;
 
%% process the high res mesh first
h_original = h;
h = 0.2;

meshname = sprintf('mesh_data%cbunny4.1',fs);

if exist([meshname '.mat'], 'file') ~= 2
    error('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end

[~, ind] = sortrows(nodeM,3);

% put the nodes in order (in z)
nodeM = nodeM(ind,:);

% create the rank for the original z
[~,R] = sort(ind);
% the new element list
elem = R(elem(:,:));

% nodeM = loaded.obj.nodeM;
% elem = loaded.obj.elem;

% fidn all the points around the same height with the lowest point
bottom_points = find(abs(nodeM(:,3)-nodeM(1,3)) < (0.05/100));

nFixed = length(bottom_points);

N = size(nodeM,1); % number of nodes

% N = size(nodeM,1); % number of nodes
indAll = 1:N;
indRemove = indAll(bottom_points(:));
indLogical = true(3,N);
indLogical(:,indRemove) = false;
indLogical = indLogical(:);
%% find elements close to constraints and make them softer
soft_points = find(abs(nodeM(:,3)-min(nodeM(:,3))) < (0.25/100));

soft_points_height_threshold = min(soft_points);

Y_list = Y * ones(size(elem,1),1);
P_list = P * ones(size(elem,1),1);
rho_list = rho * ones(size(elem,1),1);

for i = 1:size(elem,1)
    if all(elem(i,:) > soft_points_height_threshold)
        Y_list(i) = Y_list(i) / 10;
    end
end

obj = elasticTetObj(nodeM, elem);
switch material_type
    case 'linear'
        obj.SetMaterial( Y_list, P_list, rho_list, 2 * ones(size(elem,1),1) , a, b); % set the tri to linear
        %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
    case 'neo-hookean'
        obj.SetMaterial( Y_list, P_list, rho_list, 1 * ones(size(elem,1),1), a, b); % set the tri to neo-hookean
    case 'stvk'
        obj.SetMaterial( Y_list, P_list, rho_list, 3 * ones(size(elem,1),1), a, b); % set the tri to stvk
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


M = obj.M;
K = obj.StiffnessMatrix;
M = M(indLogical,indLogical);
K = K(indLogical,indLogical);
%         K = K(~ind_fix,~ind_fix); % extract the non-fixed part

[V,D] = eigs(K,M,eig_modes,'sm');
[low_eig, permutation_indices] = sort(diag(D));
V = V(:,permutation_indices);
% for i = 1:eig_modes
%     Dx(indLogical) = V(:,i)/10000;
%     obj.SetCurrentState(Dx);
%     obj.simple_vis;
%     cla
% end
high_res_low_eig = low_eig;
% v = zeros(length(Dx),1);

%%
meshname = sprintf('mesh_data%cvcoarsebunny.1',fs);

if exist([meshname '.mat'], 'file') ~= 2
    error('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end

[~, ind] = sortrows(nodeM,3);

% put the nodes in order (in z)
nodeM = nodeM(ind,:);

% create the rank for the original z
[~,R] = sort(ind);
% the new element list
elem = R(elem(:,:));

% nodeM = loaded.obj.nodeM;
% elem = loaded.obj.elem;

% fidn all the points around the same height with the lowest point
bottom_points = find(abs(nodeM(:,3)-nodeM(1,3)) < (0.05/100));

nFixed = length(bottom_points);

N = size(nodeM,1); % number of nodes

% N = size(nodeM,1); % number of nodes
indAll = 1:N;
indRemove = indAll(bottom_points(:));
indLogical = true(3,N);
indLogical(:,indRemove) = false;
indLogical = indLogical(:);
%% find elements close to constraints and make them softer
soft_points = find(abs(nodeM(:,3)-min(nodeM(:,3))) < (0.25/100));

soft_points_height_threshold = min(soft_points);

Y_list = Y * ones(size(elem,1),1);
P_list = P * ones(size(elem,1),1);
rho_list = rho * ones(size(elem,1),1);

for i = 1:size(elem,1)
    if all(elem(i,:) > soft_points_height_threshold)
        Y_list(i) = Y_list(i) / 10;
    end
end

obj = elasticTetObj(nodeM, elem);
switch material_type
    case 'linear'
        obj.SetMaterial( Y_list, P_list, rho_list, 2 * ones(size(elem,1),1) , a, b); % set the tri to linear
        %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
    case 'neo-hookean'
        obj.SetMaterial( Y_list, P_list, rho_list, 1 * ones(size(elem,1),1), a, b); % set the tri to neo-hookean
    case 'stvk'
        obj.SetMaterial( Y_list, P_list, rho_list, 3 * ones(size(elem,1),1), a, b); % set the tri to stvk
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

obj.high_res_eig_to_match = high_res_low_eig;

M = obj.M;
K = obj.StiffnessMatrix;
M = M(indLogical,indLogical);
K = K(indLogical,indLogical);


[V,D] = eigs(K,M,eig_modes,'sm');
[low_eig, permutation_indices] = sort(diag(D));
V = V(:,permutation_indices);

% if h ~= 0.5
%     V = -V;
% end
obj.eig_ratios = high_res_low_eig./low_eig;
obj.eig_targets = high_res_low_eig;
obj.eig_modes = eig_modes;
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

dirname = sprintf('sim_data%c%s_%s%d', fs, [mfilename],material_type,eig_modes);

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

        disp('simulating frame')
        disp(ti)
                u = SemiImplicit_rescaled(dt, u, obj, ~indLogical);
        
%         if (ti < 51) || (ti > 150)
        if mod(ti,40) < 20
            disp('moving up')
            positions = u(1:end/2)+obj.X;
            for bottom_i = bottom_points'
                temp = positions((3*(bottom_i-1)+1):3*bottom_i) + [0 0 1]' * 0.005 * dt;
                positions((3*(bottom_i-1)+1):3*bottom_i) = temp(1:3);
            end
            u(1:end/2) = positions-obj.X;
        else 
            disp('moving down')
            positions = u(1:end/2)+obj.X;
            for bottom_i = bottom_points'
                temp = positions((3*(bottom_i-1)+1):3*bottom_i) - [0 0 1]' * 0.005 * dt;
                positions((3*(bottom_i-1)+1):3*bottom_i) = temp(1:3);
            end
            u(1:end/2) = positions-obj.X;

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
                campos([-0.0975   -0.1267    0.1142]);
                camtarget([    0.0003    0.0006    0.0018])
                camva(6.9295);
                scal = 2000;
                offc = 0.13;
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
    
%     % plot the energy
%     for ti = 1:tsteps
%         u = trajectory(:,ti);
%         x = u(1:end/2);
%         v = u(end/2+1:end);
%         gravitational_potential(ti) = -(x+obj.X)' * obj.M * obj.externalGravity;
%         kinetic_energy(ti) = 1/2 * v' * obj.M * v;
%     end
%     cla reset;
%     plot(1:tsteps,elastic_energy,1:tsteps,kinetic_energy,1:tsteps,gravitational_potential,...
%         1:tsteps,elastic_energy+kinetic_energy+gravitational_potential)
%     print([simdir fs 'energy'],'-dpng')
    save([simdir fs 'trajectory.mat']);
else
    disp([material_type simdir fs 'trajectory.mat'])
    disp('already exist. delete it if you want to re-run')
end

name = [simdir fs 'trajectory.mat'];
close(vid);
end