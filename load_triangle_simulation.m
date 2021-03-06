function load_triangle_simulation

clear all
close all

video_name = 'comp_rect.m4v';

filename_list = {'triangle maxA1e-01_linear_longtimeTrajectory.mat','triangle maxA1e-01_linear_subTrajectory.mat','triangle maxA1e-02_linear_longtimeTrajectory.mat','triangle maxA1e-03_linear_longtimeTrajectory.mat'};
% filename_list = {'triangle maxA1e-01_linear_longtimeTrajectory.mat','triangle maxA1e-01_linear_subTrajectory.mat'};
filename_list = {'simulation_rect_maxA_5e-04.mat','simulation_rect_maxA_1e-02.mat'};
edge_colors = cell(size(filename_list));
edge_colors{1} = 'b';
edge_colors{2} = 'r';
% edge_colors{3} = 'k';
% edge_colors{4} = 'g';


fs = filesep;
trajectories = cell(size(filename_list));
positionsMs = cell(size(filename_list));
indLogicals = cell(size(filename_list));
elems = cell(size(filename_list)); % faces
for i_filename = 1:length(filename_list)
    filename = filename_list{i_filename};
    
    load(filename);
    
    trajectories{i_filename} = trajectory;
    positionsMs{i_filename} = positionsM;
    indLogicals{i_filename} = indLogical;
    elems{i_filename} = elem;
end

boundaries = cell(size(filename_list));
for i_filename = 1:length(trajectories)
    u = trajectories{i_filename}(:,1);
    
    node = u(1:end/2) + positionsMs{i_filename}(indLogicals{i_filename});
    node = transpose(reshape(node,2,[]));
    TR = triangulation(elems{i_filename}, node(:,1), node(:,2));
    boundaries{i_filename} = freeBoundary(TR)';
end

    vid = VideoWriter(video_name,'MPEG-4');
    vid.FrameRate = 60;
    open(vid);


% rate to draw the scene
sim_rate = round(1/dt);
draw_rate = round(sim_rate/vid.FrameRate);
axis_box = [-0.75 1.25 -0.5 1.5];
hf = figure;

for ti = 1:tsteps
    
    
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)
            
            for i_filename = 1:length(trajectories)
                u = trajectories{i_filename}(:,ti);
                
                
                node = u(1:end/2) + positionsMs{i_filename}(indLogicals{i_filename});
                node = transpose(reshape(node,2,[]));
                trimesh(elems{i_filename},node(:,1),node(:,2),i_filename*ones(size(node,1),1),'FaceAlpha',0.0,'EdgeAlpha',0.5,'EdgeColor',edge_colors{i_filename});
%                 bi = boundary(node(:,1),node(:,2),1); % boundary indices
%                 plot(node(boundaries{i_filename},1),node(boundaries{i_filename},2),edge_colors{i_filename});
%                 axis(axis_box)
%                 plot(node(boundaries{i_filename},1),node(boundaries{i_filename},2));
%               

                axis equal
                axis(axis_box)
                drawnow
                hold on
                legend('coarse','coarse with ghost','refined1','refined2')
                                handles_list{i_filename}.XData = node(:,1);
                                handles_list{i_filename}.YData = node(:,2);
                
            end

                frame = getframe(hf);
                writeVideo(vid,frame);

            
        end
        

    
    
    hold off
    
    
end