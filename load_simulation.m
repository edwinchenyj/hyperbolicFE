function load_simulation

clear all
close all

video_name = 'rect_comp_P_4e-01.m4v';

% filename_list = {'triangle maxA1e-01_linear_longtimeTrajectory.mat','triangle maxA1e-01_linear_subTrajectory.mat','triangle maxA1e-02_linear_longtimeTrajectory.mat','triangle maxA1e-03_linear_longtimeTrajectory.mat'};
% filename_list = {'triangle maxA1e-01_linear_longtimeTrajectory.mat','triangle maxA1e-01_linear_subTrajectory.mat'};
% filename_list = {'simulation_rect_maxA_5e-04.mat','simulation_rect_maxA_1e-02.mat','simulation_rect_maxA_1e-01.mat'};
filename_list = {'simulation_rect_maxA_1e-01_P_4e-01.mat','simulation_rect_maxA_1e-02_P_4e-01.mat'};
edge_colors = cell(size(filename_list));
edge_colors{1} = 'b';
edge_colors{2} = 'r';
% edge_colors{1} = 'b';
% edge_colors{4} = 'g';


fs = filesep;
trajectories = cell(size(filename_list));
positionsMs = cell(size(filename_list));
indLogicals = cell(size(filename_list));
elems = cell(size(filename_list)); % faces

lowests= cell(size(filename_list));
for i_filename = 1:length(filename_list)
    filename = filename_list{i_filename};
    
    S = load(filename);
    
    trajectories{i_filename} = S.trajectory;
    positionsMs{i_filename} = S.positionsM;
    indLogicals{i_filename} = S.indLogical;
    elems{i_filename} = S.elem;
    positionss{i_filename} = S.positions;
    axis_boxes{i_filename} = S.axis_box;
    dt = S.dt; % TODO: make sure all the dt are the same across files for fair comparison
    tsteps = S.tsteps;
    lowests{i_filename} = min(min(S.trajectory(1:2:end/2,:)));
end

boundaries = cell(size(filename_list));
for i_filename = 1:length(trajectories)
    u = trajectories{i_filename}(:,1);
    indLogical = indLogicals{i_filename};
    positions = positionss{i_filename};
%     node = u(1:end/2) + positionsMs{i_filename}(indLogicals{i_filename});
%     node = transpose(reshape(node,2,[]));
    positions(indLogical) = u(1:end/2) + positionsMs{i_filename}(indLogicals{i_filename});
    node = transpose(reshape(positions,2,[]));
            
    TR = triangulation(elems{i_filename}, node(:,1), node(:,2));
    boundaries{i_filename} = freeBoundary(TR)';
end

    vid = VideoWriter(video_name,'MPEG-4');
    vid.FrameRate = 60;
    open(vid);


% rate to draw the scene
sim_rate = round(1/dt);
draw_rate = round(sim_rate/vid.FrameRate);
axis_box = axis_boxes{1};
axis_box(4) = 1.5;
hf = figure;

for ti = 1:tsteps
    
    
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)
            
            for i_filename = 1:length(trajectories)
                u = trajectories{i_filename}(:,ti);
                
                indLogical = indLogicals{i_filename};
                positions = positionss{i_filename};
                positionsM = positionsMs{i_filename};
                positions(indLogical) = u(1:end/2) + positionsM(indLogical);
                node = transpose(reshape(positions,2,[]));
            
                
%                 node = u(1:end/2) + positionsMs{i_filename}(indLogicals{i_filename});
%                 node = transpose(reshape(node,2,[]));
                trimesh(elems{i_filename},node(:,2),node(:,1),i_filename*ones(size(node,1),1),'FaceAlpha',0.0,'EdgeAlpha',0.5,'EdgeColor',edge_colors{i_filename});
%                 bi = boundary(node(:,1),node(:,2),1); % boundary indices
%                 plot(node(boundaries{i_filename},2),node(boundaries{i_filename},1),edge_colors{i_filename});
%                 axis(axis_box)
%                 plot(node(boundaries{i_filename},2),node(boundaries{i_filename},1));
%               

                axis equal
                axis(axis_box)
                hl = refline(0,lowests{i_filename});
                hl.Color = edge_colors{i_filename};
                drawnow
                hold on
%                 legend('coarse','coarse with ghost','refined1','refined2')
%                                 handles_list{i_filename}.XData = node(:,1);
%                                 handles_list{i_filename}.YData = node(:,2);
                
            end

                frame = getframe(gca);
                writeVideo(vid,frame);

            
        end
        

    
    
    hold off
    
    
end