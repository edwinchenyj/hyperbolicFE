function load_triangle_simulation

clear all
close all

draw = true;
rerun_flag = true;
save_state = true;
draw = true;

fs = filesep;
mass_flag = true;

type = 'triangle';

axis_box = [-1 1.5 -0.5 1];
% options = optimoptions('fsolve','Algorithm','levenberg-marquardt');
options = optimoptions('fsolve','TolFun',1.e-9,'TolX',1.e-9,'Display','final');

maxA_list = [0.01 0.001]';
%         meshname = sprintf('simData%ctrimeshes%c%s maxA%.e',fs,fs,type, maxA);
trajectories = cell(size(maxA_list));
positionsMs = cell(size(maxA_list));
indLogicals = cell(size(maxA_list));
elems = cell(size(maxA_list)); % faces
% meshname = sprintf('simData%ctrimeshes%c%s maxA%.e tritool%s',fs,fs,type, maxA);
for i_maxA = 1:length(maxA_list)
    maxA = maxA_list(i_maxA);
    meshname = sprintf('simData%ctrimeshes%c%s maxA%.e%s',fs,fs,type, maxA,'triangle');
    if exist([meshname '.mat'], 'file') ~= 2
        disp('mesh does not exist')
        
    else
        load([meshname '.mat'], 'nodeM', 'elem');
        
    end
    
    
    elem(:,[1 3]) = elem(:,[3 1]);
    
    N = size(nodeM,1);
    
    Y = 10; % Young's modululs
    P = 0.48; % Poisson ratio
    rho = 1; % density
    a = 0.0; % rayleigh damping
    b = 0.0;
    
    filename = sprintf('simData%c%s maxA%.e', fs, type, maxA);
    filename = [filename '_linear_longtime'];
    
    
    fname = [filename 'Trajectory'];
    
    load(fname,'trajectory','tsteps','dt','positionsM','indLogical');
    
    trajectories{i_maxA} = trajectory;
    positionsMs{i_maxA} = positionsM;
    indLogicals{i_maxA} = indLogical;
    elems{i_maxA} = elem;
end

boundaries = cell(size(maxA_list));
for iA = 1:length(trajectories)
    u = trajectories{iA}(:,1);
    
    
    node = u(1:end/2) + positionsMs{iA}(indLogicals{iA});
    node = transpose(reshape(node,2,[]));
    TR = triangulation(elems{iA}, node(:,1), node(:,2));
    boundaries{iA} = freeBoundary(TR)';
end

if draw
    
    vidname = strcat('linear_triangle_comp_longtime','.avi');
    vid = VideoWriter(vidname);
    vid.FrameRate = 60;
    open(vid);
end

edge_colors = cell(size(maxA_list));
edge_colors{1} = 'b';
edge_colors{2} = 'r';
% rate to draw the scene
sim_rate = round(1/dt);
draw_rate = round(sim_rate/vid.FrameRate);
axis_box = [-0.75 1.25 -0.5 1.5];
hf = figure;

for ti = 1:tsteps
    
    
    if(draw)
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)
            
            for iA = 1:length(trajectories)
                u = trajectories{iA}(:,ti);
                
                
                node = u(1:end/2) + positionsMs{iA}(indLogicals{iA});
                node = transpose(reshape(node,2,[]));
%                 trimesh(elems{iA},node(:,1),node(:,2),iA*ones(size(node,1),1),'FaceAlpha',0.0,'EdgeAlpha',0.5,'EdgeColor',edge_colors{iA});
%                 bi = boundary(node(:,1),node(:,2),1); % boundary indices
                plot(node(boundaries{iA},1),node(boundaries{iA},2),edge_colors{iA});
%                 axis(axis_box)
                
                axis equal
                axis(axis_box)
                drawnow
                hold on
                %                 handles_list{iA}.XData = node(:,1);
                %                 handles_list{iA}.YData = node(:,2);
                
            end
            if save_state
                frame = getframe(hf);
                writeVideo(vid,frame);
            end
            
        end
        
    end
    
    
    hold off
    
    
end