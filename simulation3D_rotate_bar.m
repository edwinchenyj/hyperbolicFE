function name = simulation3D_rotate_bar(varargin)
close all
clf

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

if (exist([simdir fs 'trajectory.mat'], 'file') ~= 2)
  
        if exist([meshfile '.mat'], 'file') ~= 2
            error('mesh does not exist')
        else
            load([meshfile '.mat'], 'nodeM', 'elem');
            
        end
%% set fixed points and points to be rotated
% fix nFixed number of highest point

[~, ind] = sortrows(nodeM,1);

% put the nodes in order (in z)
nodeM = nodeM(ind,:);


left_points = find(abs(nodeM(:,1)-min(nodeM(:,1))) < 0.01);
right_points = find(abs(nodeM(:,1)-max(nodeM(:,1))) < 0.01);

% create the rank for the original z
[~,R] = sort(ind);
% the new element list
elem = R(elem(:,:));

obj = elasticTetObj(nodeM, elem);
switch material_type
    case 'linear'
        obj.SetMaterial( Y, P, rho, 2, a, b); % set the tri to linear
        %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
    case 'neo-hookean'
        obj.SetMaterial( Y, P, rho, 1, a, b); % set the tri to neo-hookean
    case 'stvk'
        obj.SetMaterial( Y, P, rho, 3, a, b); % set the tri to neo-hookean
end

switch gravity
    case 'on'
        obj.gravity_on = true;
        obj.calculateGravity;
end
%% get logical indices for the left and right nodes

indLeft = left_points;
indRight = right_points;

nFixed = length(indLeft) + length(indRight);

N = size(nodeM,1); % number of nodes
indAll = 1:N;
indRemove = indAll([indLeft(:); indRight(:)]);
indLogical = logical(ones(3,N));
indLogical(:,indRemove) = logical(0);
indLogical = indLogical(:);

indLeftLogical = logical(zeros(3,N));
indLeftLogical(:,indLeft) = logical(1);
indLeftLogical = indLeftLogical(:);

indRightLogical = logical(zeros(3,N));
indRightLogical(:,indRight) = logical(1);
indRightLogical(:);

% Dx_free = Dx(indLogical); v_free = v(indLogical);

switch gravity
    case 'on'
        obj.gravity_on = true;
        obj.calculateGravity;
end

%
Dx = 0*rand(obj.Dim*N,1); % displacement field. set zero for rest state
obj.SetCurrentState(Dx);
% ha = obj.init_vis;
% obj.simple_vis(obj.vis_handle);
    ha = obj.init_vis;
    obj.simple_vis(obj.vis_handle);
    
    v = zeros(length(Dx),1);
    u = [Dx; v];
    fig = gcf;
    fig.Position(1:2) = [100 50];
    fig.Position(3:4) = fig.Position(3:4) *2;
    
    vidname = strcat(simdir,fs,'video.avi');
    %     vid = VideoWriter(vidname,'Uncompressed AVI');
    vid = VideoWriter(vidname);
    vid.FrameRate = 50;
    open(vid);
    % rate to draw the scene
    sim_rate = round(1/dt);
    draw_rate = round(sim_rate/vid.FrameRate);
    
    
    kinetic_energy = zeros(1,tsteps);
    
    patch_color = [0.8,0.9,0.9];
    trajectory = zeros(size(u,1),tsteps);
    sim_res = zeros(1,tsteps);
    sim_gradient_size = zeros(1,tsteps);
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
    annotation('textbox','String',sim_annotation,'FitBoxToText','on','LineStyle','none','FontSize',12);
    
    
%% define the rotation transform

angular_v = 1;

left_transform = [1 0 0 0; 0 1 0 -1/6/5*10; 0 0 1 -1/6/5*10; 0 0 0 1];
inv_left_transform = [1 0 0 0; 0 1 0 1/6/5*10; 0 0 1 1/6/5*10; 0 0 0 1];

right_transform = [1 0 0 0; 0 1 0 -1/6/5*10; 0 0 1 -1/6/5*10; 0 0 0 1];
inv_right_transform = [1 0 0 0; 0 1 0 1/6/5*10; 0 0 1 1/6/5*10; 0 0 0 1];

rotation_transform = [ 1 0 0 0; 0 cos(angular_v * dt) -sin(angular_v * dt) 0; 0 sin(angular_v * dt) cos(angular_v * dt) 0; 0 0 0 1];
inv_rotation_transform = [ 1 0 0 0; 0 cos(-angular_v * dt) -sin(-angular_v * dt) 0; 0 sin(-angular_v * dt) cos(-angular_v * dt) 0; 0 0 0 1];

patch_color = [0.8,0.9,0.5];

    for ti = 1:tsteps
        disp('frame'); disp(ti)
        trajectory(:,ti) = u;
        
        [u, res, step_size]= solver(dt, u, obj,~indLogical);
        
        if ti < 201
            % apply the rotation
            positions = u(1:end/2)+obj.X;
            for left_i = indLeft'
                temp = inv_left_transform * rotation_transform *left_transform * [positions((3*(left_i-1)+1):3*left_i); 1];
                positions((3*(left_i-1)+1):3*left_i) = temp(1:3);
            end
            
%             for right_i = indRight'
%                 temp = inv_right_transform * inv_rotation_transform *right_transform * [positions((3*(right_i-1)+1):3*right_i); 1];
%                 positions((3*(right_i-1)+1):3*right_i) = temp(1:3);
%             end
            u(1:end/2) = positions-obj.X;
        end
        
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)
            axis equal
            %                 axis auto
            cla
            
            
            obj.simple_vis(obj.vis_handle,patch_color);
            %
            ah = gca;
            ph = ah.Children;
            ph.FaceLighting = 'flat';
            shading interp
            ph.EdgeLighting = 'none';
            ph.EdgeAlpha = 0;
            whitebg('white')
            if any(strfind(meshfile,'small_bar'))
                campos([-1.1513   -1.6081    1.4854]);
                camtarget([0.0810   -0.0021    0.0680])
                camva(6.9295);
            elseif any(strfind(meshfile,'octopus'))
                campos([   -3.7639   -4.9106    3.4267]);
                camtarget([    0.0603    0.0732   -0.2002])
                camva(6.9295);
            elseif any(strfind(meshfile,'armadillo'))
                campos([2.4703  -20.8381    5.3614]);
                camtarget([0.7785   -0.6424    0.9054])
                camva(6.9295);
            elseif any(strfind(meshfile,'horse'))
                campos([    -12.6549   -9.8499   17.1693]);
                camtarget([    -0.0506   -0.0023    0.0166])
%                 camva(7.3921);
                camup([    0.5763    0.4503    0.6820]);
            end
            hold on
            colormap(bone);
            % obj.patch_handle.FaceVertexCData = ones(63,3) *0.5;
            facecolor = winter(6);
            %                 facecolor = facecolor(2,:);
            obj.patch_handle.FaceVertexCData = repmat(facecolor(5,:),N,1);
            obj.patch_handle.FaceAlpha = 0.2;
            obj.patch_handle.EdgeColor = [0.1 0.1 0.1];
            obj.patch_handle.EdgeAlpha = 0.8;
            obj.patch_handle.Marker = 'none';
            facecolor2 = spring(6);
            obj.patch_handle.MarkerFaceColor = facecolor2(4,:);
            
            fig.Children.Visible = 'off';
            fig.Children.Clipping = 'off';
            fig.Children.Projection = 'perspective';
            lighting flat
            material dull
            lh = camlight('headlight');
            lh.Style = 'local';
            %
            %             dim = [.7 .0 .7 .5];
            
            frame = getframe(fig);
            writeVideo(vid,frame);
            
            
        end
    end
    
    % plot the energy
    for ti = 1:tsteps
        u = trajectory(:,ti);
        x = u(1:end/2);
        v = u(end/2+1:end);
        kinetic_energy(ti) = 1/2 * v' * obj.M * v;
    end
    cla reset;
    Mke = max(kinetic_energy);
    plot(1:tsteps,kinetic_energy/Mke)
    print([simdir fs 'energy'],'-dpng')
    figure
    plot(sim_res);
    title('residual');
    figure
    plot(sim_gradient_size);
    title('gradient size');
    
    save([simdir fs 'trajectory.mat']);
    close(gcf)
    close(gcf)
    close(gcf)
    close(vid);
    name = [simdir fs 'trajectory.mat'];
else
    disp([material_type simdir fs 'trajectory.mat'])
    disp('already exist. delete it if you want to re-run')
    name = [simdir fs 'trajectory.mat'];
    
end
