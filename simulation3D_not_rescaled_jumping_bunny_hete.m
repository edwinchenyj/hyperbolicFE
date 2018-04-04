function name = simulation3D_not_rescaled_jumping_bunny_hete(varargin)
close all
clf

fs = filesep;

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
            
    end
    i_arg = i_arg + 1;
end

tsteps = T/dt;


%%
scriptdir = script_directory(mfilename);
mkdir(scriptdir);
simdir = [scriptdir,fs,sim_directory(varargin)];
mkdir(simdir);

if (exist([simdir fs 'trajectory.mat'], 'file') ~= 2)
    

    
    %%  original mesh
    % h = h_original;
    if strcmp(mesh_shape,'coarsebunny')
        meshname = sprintf('mesh_data%ccoarsebunny.1',fs);
    elseif strcmp(mesh_shape,'vcoarsebunny')
        meshname = sprintf('mesh_data%cvcoarsebunny.1',fs);
    end
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
bottom_points = find(abs(nodeM(:,3)-nodeM(1,3)) < (0.01/100));

nFixed = length(bottom_points);

N = size(nodeM,1); % number of nodes

% N = size(nodeM,1); % number of nodes
indAll = 1:N;
indRemove = indAll(bottom_points(:));
indLogical = true(3,N);
indLogical(:,indRemove) = false;
indLogical = indLogical(:);
%% find elements close to constraints and make them softer
soft_points = find(abs(nodeM(:,3)-min(nodeM(1,3))) < (0.4/100));

soft_points_height_threshold = max(soft_points);

Y_list = Y * ones(size(elem,1),1);
P_list = P * ones(size(elem,1),1);
rho_list = rho * ones(size(elem,1),1);

for i = 1:size(elem,1)
    if all(elem(i,:) < soft_points_height_threshold)
        Y_list(i) = Y_list(i)/10;
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
    
    ha = obj.init_vis;
    obj.simple_vis(obj.vis_handle);
    
    v = zeros(length(Dx),1);
    axis equal
    u = [Dx; v];
    fig = gcf;
    
    fig.Children.Visible = 'off';
    fig.Children.Clipping = 'off';
    fig.Children.Projection = 'perspective';
    fig.Position(1:2) = [100 50];
    fig.Position(3:4) = fig.Position(3:4) *2;
    
    
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
    
    kinetic_energy = zeros(1,tsteps);
    
    
    patch_color = [0.8,0.9,0.5];
    
    
    
    trajectory = zeros(size(u,1),tsteps);
    for ti = 1:tsteps
        trajectory(:,ti) = u;
        %         elastic_energy(ti) = obj.ElasticEnergy;
        
        switch solver
            case 'SI'
                u = SemiImplicit(dt, u, obj,~indLogical);
            case 'ERE'
                u = ERE(dt, u, obj,~indLogical);
            case 'BE'
                u = BackwardEuler_line_search(dt, u, obj,~indLogical);
%                 u = eig_BackwardEuler(dt, u, obj, ~indLogical);
            case 'IM'
                u = ImplicitMid_line_search(dt, u, obj,~indLogical);
            otherwise
                error('no solver')
        end
        
%         if (ti < 51) || (ti > 150)
        if mod(ti,40) < 20
            disp('moving up')
            positions = u(1:end/2)+obj.X;
            for bottom_i = bottom_points'
                temp = positions((3*(bottom_i-1)+1):3*bottom_i) + [0 0 1]' * 0.005 * dt;
                positions((3*(bottom_i-1)+1):3*bottom_i) = temp(1:3);
            end
            u(1:end/2) = positions-obj.X;
%         elseif ti < 40
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
            campos([-0.0975   -0.1267    0.1142]);
            camtarget([    0.0003    0.0006    0.0018])
                camva(6.9295);
                hold on
                colormap(bone);
                % obj.patch_handle.FaceVertexCData = ones(63,3) *0.5;
                facecolor = winter(6);
                %                 facecolor = facecolor(2,:);
                obj.patch_handle.FaceVertexCData = repmat(facecolor(5,:),N,1);
                obj.patch_handle.FaceAlpha = 0.2;
                obj.patch_handle.EdgeColor = [0.1 0.1 0.1];
                obj.patch_handle.EdgeAlpha = 0.8;
%                 obj.patch_handle.Marker = 'none';
%                 facecolor2 = spring(6);
%                 obj.patch_handle.MarkerFaceColor = facecolor2(4,:);
                
                fig.Children.Visible = 'off';
                fig.Children.Clipping = 'off';
                fig.Children.Projection = 'perspective';
                lighting flat
                material dull
                lh = camlight('headlight');
                lh.Style = 'local';
                
                
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
        kinetic_energy(ti) = 1/2 * v' * obj.M * v;
    end
    cla reset;
    Mke = max(kinetic_energy);
    plot(1:tsteps,kinetic_energy/Mke)
    print([simdir fs 'energy'],'-dpng')
    save([simdir fs 'trajectory.mat']);
    if draw
    close(vid);
    end
    name = [simdir fs 'trajectory.mat'];
else
    disp([material_type simdir fs 'trajectory.mat'])
    disp('already exist. delete it if you want to re-run')
    name = [simdir fs 'trajectory.mat'];

end


end