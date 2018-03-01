function filename_list=releasing_bar_video_3comparsion_polyfit_script

clear all
close all

video_name = 'CG_releasing_3bar_Y1000_neo-hookean_DAC.avi';

fs = filesep;

% compare 2 simulations
filename_list = cell(1,3);

%% first comparision
draw = true;
rerun_flag = true;
save_state = true;
dt = [1/100];
T = [2];
mesh_shape = 'small_bar';
simulation_type = 'CG';
solver = 'ERE';

Y = [1e3];
P = [0.45];
rho = [1000];
a = [0.01];
b = [0.01];
constraint = 'n';
material_type = 'neo-hookean';
gravity = 'off';
% DGeta = 1e6;
axis_box = [-1 1 -1 1 -1 1]/4;
deformation_scale = 'n';


xlim = axis_box(1:2);
ylim = axis_box(3:4);
zlim = axis_box(5:6);

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

h = [0.5];


dirname = sprintf('sim_data%csimulation3D_releasing_bar_%s_%s', fs,material_type, mesh_shape);
if isDG
    simdir = strcat(dirname,fs,solver,'_',...
        simulation_type,...
        '_constraint_',constraint,...
        '_h',num2str(h),...
        '_Y',num2str(Y),...
        '_P',num2str(P),...
        '_rho',num2str(rho),...
        '_a',num2str(a),...
        '_b',num2str(b),...
        '_dt',num2str(dt),...
        '_eta',num2str(DGeta),...
        '_def-scl',num2str(deformation_scale));
else
  
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
    
end
filename_list{1} = [simdir fs 'trajectory.mat'];


%% second comparison

simulation_type = 'CG';
solver = 'ERE';
h = [0.5];
polyfit_modes = 2;
is_polyfit = true;
dirname = sprintf('sim_data%csimulation3D_releasing_bar_DAC_%s_%s', fs,material_type, mesh_shape);


if isDG
    simdir = strcat(dirname,fs,solver,'_',...
        '_polyfit_',...
        simulation_type,...
        '_constraint_',constraint,...
        '_h',num2str(h),...
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
        '_constraint_',constraint,...
        '_h',num2str(h),...
        '_Y',num2str(Y),...
        '_P',num2str(P),...
        '_rho',num2str(rho),...
        '_a',num2str(a),...
        '_b',num2str(b),...
        '_dt',num2str(dt),...
        '_def-scl',num2str(deformation_scale));
    
end

filename_list{2} = [simdir fs 'trajectory.mat'];


%% second comparison

simulation_type = 'CG';
solver = 'ERE';
h = [0.2];
polyfit_modes = 2;
is_polyfit = true;
dirname = sprintf('sim_data%csimulation3D_releasing_bar_%s_%s', fs,material_type, mesh_shape);


if isDG
    simdir = strcat(dirname,fs,solver,'_',...
        '_polyfit_',...
        simulation_type,...
        '_constraint_',constraint,...
        '_h',num2str(h),...
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
        '_constraint_',constraint,...
        '_h',num2str(h),...
        '_Y',num2str(Y),...
        '_P',num2str(P),...
        '_rho',num2str(rho),...
        '_a',num2str(a),...
        '_b',num2str(b),...
        '_dt',num2str(dt),...
        '_def-scl',num2str(deformation_scale));
    
end

filename_list{3} = [simdir fs 'trajectory.mat'];


%%

edge_colors = num2cell(hsv(4),2); % create a list of different colours
trajectories = cell(size(filename_list));
positionsMs = cell(size(filename_list));
indLogicals = cell(size(filename_list));
elems = cell(size(filename_list)); % faces
positionss = cell(size(filename_list));
axis_boxes = cell(size(filename_list));
lowests= cell(size(filename_list));
h = cell(size(filename_list));

vid = VideoWriter(video_name);
vid.FrameRate = 50;
open(vid);



% visualizing for CG
if strcmp(simulation_type(1:2),'CG')
    
    for i_filename = 1:length(filename_list)
        filename = filename_list{i_filename};
        
        S = load(filename);
        
        trajectories{i_filename} = S.trajectory;
        positionsMs{i_filename} = S.obj.X;
        indLogicals{i_filename} = S.indLogical;
        elems{i_filename} = S.obj.elem;
        positionss{i_filename} = S.obj.X;
        axis_boxes{i_filename} = S.axis_box;
        Ns(i_filename) = S.N;
        dt = S.dt; % TODO: make sure all the dt are the same across files for fair comparison
        tsteps = S.tsteps;
        lowests{i_filename} = min(min(S.trajectory(1:2:end/2,:)));
    end
    
    % boundary for the mesh in CG sim
    boundaries = cell(size(filename_list));
    for i_filename = 1:length(trajectories)
        u = trajectories{i_filename}(:,1);
        indLogical = indLogicals{i_filename};
        positions = positionss{i_filename};
        
        positions = u(1:end/2) + positionsMs{i_filename};
        node = transpose(reshape(positions,3,[]));
        
%         TR = triangulation(elems{i_filename}, node(:,1), node(:,2));
%         boundaries{i_filename} = freeBoundary(TR)';
    end
    
    
    % rate to draw the scene
    sim_rate = round(1/dt);
    draw_rate = round(sim_rate/vid.FrameRate);
%     axis_box = axis_boxes{1};
%     axis_box(4) = 1;
%     axis_box(1) = -0.5;
%     axis_box(2) = 1.5;
    
%     hf = figure;
clf
ha = S.obj.init_vis;
S.obj.simple_vis(S.obj.vis_handle);
    fig = gcf;

fig.Children.Visible = 'off';
fig.Children.Clipping = 'off';
fig.Children.Projection = 'perspective';
fig.Position(3:4) = fig.Position(3:4) *1;
fig.Position(1:2) = [0,0];
                ccmap = gray - 0.8;
                ind = ccmap < 0;
                ccmap(ind) = 0;
skymap = colormap(ccmap); %We'll use the colormap "winter" which is defined by matlab.

% Here's a gradient I made up that's a bit more horizon-like. You can
% experiment.
% flipud([logspace(-0.5,0,1000)',logspace(-0.5,0,1000)',logspace(-0.001,0,1000)']);

[skyX,skyY,skyZ] = sphere(200); %create the surface data for a sphere.
cla
hold on
%                     lh = camlight('headlight');
%                 lh.Style = 'local';

                
    for ti = 1:tsteps
        
        
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)
            hold on
            scal = 5;
            offc = 0.05;
            sky = surf(scal*(skyX-offc),scal*(skyY-offc),scal*(skyZ-offc),'LineStyle','none','FaceColor','flat','FaceLighting','gouraud'); %Make a sphere object.
            
            for i_filename = 1:length(trajectories)

                u = trajectories{i_filename}(:,ti);
                
                indLogical = indLogicals{i_filename};
                positions = positionss{i_filename};
                positionsM = positionsMs{i_filename};
                positions = u(1:end/2) + positionsM;
                switch i_filename
                    case 1
                        positions(2:3:end) = positions(2:3:end) - 0.08;
                        positions(3:3:end) = positions(3:3:end) + 0.04;
                    case 2
                        positions(2:3:end) = positions(2:3:end) + 0.07;
                        positions(3:3:end) = positions(3:3:end) + 0.04;
                    case 3
                        positions(2:3:end) = positions(2:3:end) + 0.22;
                        positions(3:3:end) = positions(3:3:end) + 0.04;
                end
                node = transpose(reshape(positions,3,[]));
                tri1 = surftri(node,elems{i_filename});
                h{i_filename} = trimesh(tri1,node(:,1),node(:,2),node(:,3));
%                 triplot(elems{i_filename},node(:,1),node(:,2),'Color',edge_colors{i_filename});
                
                axis equal
                axis(axis_box)
                % draw a horizontal reference line for droop test
                %                 hl = refline(0,lowests{i_filename});
                %                 hl.Color = edge_colors{i_filename};
                
%                 ah = gca;
%                 phlist = findobj(ah,'Type','patch');
%                 ph = phlist(i_filename);
                h{i_filename}.FaceLighting = 'flat';
                h{i_filename}.FaceAlpha = 0.3;
                h{i_filename}.FaceColor = 'interp';
                %                 ph.EdgeLighting = 'gouraud';
                %                 ph.BackFaceLighting = 'lit';

                %                 ph.EdgeColor = [0.9 0.9 1];

                %                 lh2 = camlight('headlight');
                %                 lh.Position(2) = 3;
                %                 lh.Position(3) = 3;
                cmap = [0.8 0.9 1];
                h{i_filename}.FaceVertexCData = repmat(edge_colors{i_filename},Ns(i_filename),1);
                %                 colormap(autumn);
                
                shading faceted
                h{i_filename}.EdgeColor = edge_colors{i_filename};
%                 ph.EdgeLighting = 'none';
                h{i_filename}.EdgeAlpha = 0.7;
                %                 ph.LineWidth = 0.25;
                %                 ph.AlignVertexCenters = 'on';
%                 whitebg(cmap)
                
                %                 origPos = [0,0,0]';
                % %                 campos(origPos);
                %                 origTarget = [0,0,-10000];
                %                 camtarget(origTarget);
                %                 camva(40)
                ha.Visible = 'off';
                ha.Clipping = 'off';
                ha.Projection = 'perspective';
                campos([-1.1513   -1.6081    0.9854]);
                camtarget([0.0810   -0.0021    0.0680])
                camva(6.9295);
                lighting gouraud
                material metal

                colormap(ccmap)
                hold on
%                 add_shadow()

                drawnow
            end
            
            lh = camlight('headlight');
            lh.Style = 'local';
            lh.Position(3) = lh.Position(3) * 2;
            ground = [];
            color = [0.21 0.21 0.21];
            nudge = 0.12;
            background_color = [1 1 1];
            fade = 'local';
            T = [h{:}]';
            if isempty(ground)
                minZ = inf;
                for t = T'
                    V = t.Vertices;
                    minZ = min([V(:,3);minZ]);
                end
                ground = [0 0 -1 -nudge];
            end
            for i_filename = 1:length(trajectories)
                t = h{i_filename};
                t.EdgeColor = edge_colors{i_filename};
                t.FaceAlpha = 0.5;
                t.EdgeAlpha = 0.6;
            
            end
            h_shadow = [];
            M = zeros([0,0,0]);
            L = lh;
            for t = T'
                V = t.Vertices;
                for l = L'
                    % plane equation
                    % 0 = ax + by + cz + d
                    light_pos = [l.Position strcmp(l.Style,'local')];
                    d = ground * light_pos';
                    shadow_mat = d*eye(4) - light_pos'*ground;
                    U = [V ones(size(V,1),1)]*shadow_mat';
                    U = bsxfun(@rdivide,U(:,1:3),U(:,4));
                    
                    hold on;
                    tsh = trisurf(t.Faces,U(:,1),U(:,2),U(:,3), ...
                        'FaceColor',color, ...
                        'DiffuseStrength',0,'SpecularStrength',0, ...
                        'AmbientStrength',1, ...
                        'EdgeColor','none');
                    hold off;
                    switch fade
                        case {'local','infinite'}
                            D = matrixnormalize( ...
                                sum(bsxfun(@times,tsh.Vertices(:,1:2),l.Position(1:2)),2));
                            switch fade
                                case 'infinite'
                                    D = 1.0-D;
                            end
                            tsh.FaceVertexCData = ...
                                bsxfun(@plus,color,bsxfun(@times,D,background_color-color)/3-0.2);
                            tsh.FaceColor = 'interp';
                    end
                    h_shadow = [h_shadow;tsh];
                    M(:,:,end+1) = shadow_mat;
                end
            end
%             add_shadow();
%             legend('fine','coarse with 1 polyfit','coarse with 2 polyfit','coarse','Location','southeast')
            frame = getframe(gca);
            writeVideo(vid,frame);
        end
        
        cla
        hold off
        
    end
    
    % visualizing for DG
elseif(strcmp(simulation_type(1:2),'DG'))
    
    obj = cell(size(filename_list));
    
    for i_filename = 1:length(filename_list)
        filename = filename_list{i_filename};
        
        S = load(filename);
        
        trajectories{i_filename} = S.trajectory;
        positionsMs{i_filename} = S.positionsM;
        indLogicals{i_filename} = S.indLogical;
        %         elems{i_filename} = S.DGelem;
        positionss{i_filename} = S.positions;
        axis_boxes{i_filename} = S.axis_box;
        dt = S.dt; % TODO: make sure all the dt are the same across files for fair comparison
        tsteps = S.tsteps;
        lowests{i_filename} = min(min(S.trajectory(1:2:end/2,:)));
        obj{i_filename} = S.obj;
        
    end
    
    % rate to draw the scene
    sim_rate = round(1/dt);
    draw_rate = round(sim_rate/vid.FrameRate);
    axis_box = axis_boxes{1};
    axis_box(4) = 1.5;
    hf = figure;
    ha = axes;
    axis(axis_box)
    axis equal
    xlim = ha.XLim;
    ylim = ha.YLim;
    for i_filename = 1:length(filename_list)
        obj{i_filename}.init_vis(ha);
    end
    
    for ti = 1:tsteps
        
        
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)
            
            for i_filename = 1:length(trajectories)
                u = trajectories{i_filename}(:,ti);
                
                indLogical = indLogicals{i_filename};
                positions = positionss{i_filename};
                positionsM = positionsMs{i_filename};
                positions(indLogical) = u(1:end/2) + positionsM(indLogical);

                obj{i_filename}.SetCurrentDGState(positions - positionsM);
                obj{i_filename}.DG_current_vis(obj{i_filename}.vis_handle,edge_colors{i_filename});
                % draw a horizontal reference line for droop test
                %                 hl = refline(0,lowests{i_filename});
                %                 hl.Color = edge_colors{i_filename};
                
                obj{i_filename}.vis_handle.XLim = xlim;
                obj{i_filename}.vis_handle.YLim = ylim;
                
                drawnow
                hold on
            end
            
            frame = getframe(gca);
            writeVideo(vid,frame);
            
            
        end
        
        hold off
        
    end
end