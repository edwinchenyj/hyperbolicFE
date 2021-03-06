function filename_list=bar_video_4comparsion_reduced_rescaled_script

clear all
close all

video_name = 'CG_bar_Y10000_neo-hookean_reduced_ERE_rescaled_coarse.avi';

fs = filesep;

axis_box = [-1 1 -1 1 -1 1]/4;
% compare 2 simulations
filename_list = cell(1,4);

%% first comparision
filename_list{1} = 'sim_data\simulation3D_reduced_ERE_rescaled_bar6_neo-hookean_small_bar\ERE_CG_constraint_n_h0.5_Y10000_P0.45_rho1000_a0_b0.01_dt0.01_def-scln_redu20\trajectory.mat';

%% first comparision

filename_list{2} = 'sim_data\simulation3D_reduced_ERE_bar6_neo-hookean_small_bar\ERE_CG_constraint_n_h0.5_Y10000_P0.45_rho1000_a0_b0.01_dt0.01_def-scln_redu20\trajectory.mat';


%% second comparison

filename_list{3} = 'sim_data\simulation3D_not_reduced_ERE_bar6_neo-hookean_small_bar\ERE_CG_constraint_n_h0.5_Y10000_P0.45_rho1000_a0_b0.01_dt0.01_def-scln\trajectory.mat';


%% third comparison
filename_list{4} = 'sim_data\simulation3D_not_reduced_ERE_bar6_neo-hookean_small_bar\ERE_CG_constraint_n_h0.2_Y10000_P0.45_rho1000_a0_b0.01_dt0.01_def-scln\trajectory.mat';



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

              tsteps = 200;  
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
                        positions(2:3:end) = positions(2:3:end) - 0.3;
                        positions(3:3:end) = positions(3:3:end) + 0.04;
                    case 2
                        positions(2:3:end) = positions(2:3:end) - 0.2;
                        positions(3:3:end) = positions(3:3:end) + 0.04;
                    case 3
                        positions(2:3:end) = positions(2:3:end) - 0.1;
                        positions(3:3:end) = positions(3:3:end) + 0.04;
                    case 4
                        positions(2:3:end) = positions(2:3:end) + 0.0;
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
                campos([-1.3475    3.5906    2.0109]);
                camtarget([-0.0167   -0.0655    0.1104])
%                 campos([0    2.0   2]);
%                 camtarget([0   -0.15    0.1])
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
            frame = getframe(gcf);
            writeVideo(vid,frame);
        end
        
        cla
        hold off
        
    end
    
end
    
