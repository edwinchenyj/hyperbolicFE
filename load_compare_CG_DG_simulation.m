function load_compare_CG_DG_simulation


close all

video_name = 'CG_vs_DGBZ_eta_1e-1_triangle_maxA1e-1.avi';

fs = filesep;

draw = true;
rerun_flag = true;
save_state = true;
draw = true;

dt = 1/120;
tsteps = 120*1;

fs = filesep;

mesh_shape = 'triangle';
maxA = 1e-1;
% simulation_type = 'DGBZ';
%
% % set the DG flag base on simulation type
% switch simulation_type(1:2)
%     case 'DG'
%         isDG = true;
%     otherwise
%         isDG = false;
% end

DGeta = 1e-1;
solver = 'SI';
constraints = 1; % types of constraint
% 1: free

deformation_mode_number = 1;
switch maxA
    case 0.1
        deformation_scale_factor = 2;
    case 0.01
        deformation_scale_factor = -2; % there is a sign change when maxA = 0.01, 0.001
    case 0.001
        deformation_scale_factor = -2; % there is a sign change when maxA = 0.01, 0.001
end
Y = 100; % Young's modululs
P = 0.48; % Poisson ratio
rho = 1; % density
a = 0.0; % rayleigh damping
b = 0.00;
material = 'linear'; % choice: 'linear', 'neo-hookean'

axis_box = [-1 1.5 -0.5 1];

simulation_type_list = {'CG','DGBZ'};
filename_list = cell(1,length(simulation_type_list));

for i_simulation_type = 1:length(simulation_type_list)
    simulation_type = simulation_type_list(i_simulation_type);
    dirname = sprintf('sim_data%c%s_%s_maxA_%.d', fs, material, mesh_shape, maxA);
    if ~strcmp(simulation_type,'CG')
        simdir = strcat(dirname,fs,solver,'_',simulation_type,'_Y',num2str(Y),'_P',num2str(P),'_dt',num2str(dt),'_eta',num2str(DGeta));
    else
        simdir = strcat(dirname,fs,solver,'_',simulation_type,'_Y',num2str(Y),'_P',num2str(P),'_dt',num2str(dt));
    end
    filename_list(i_simulation_type) = strcat(simdir, fs, 'trajectory.mat');
end
edge_colors = num2cell(hsv(length(simulation_type_list)),2); % create a list of different colours

trajectories = cell(size(filename_list));
positionsMs = cell(size(filename_list));
indLogicals = cell(size(filename_list));
elems = cell(size(filename_list)); % faces
positionss = cell(size(filename_list));
axis_boxes = cell(size(filename_list));
obj = cell(size(filename_list));
lowests= cell(size(filename_list));


vid = VideoWriter(video_name,'MPEG-4');
vid.FrameRate = 60;
open(vid);



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

% boundary for the mesh in CG sim
boundaries = cell(size(filename_list));
for i_filename = 1:length(trajectories)
    u = trajectories{i_filename}(:,1);
    indLogical = indLogicals{i_filename};
    positions = positionss{i_filename};
    
    positions(indLogical) = u(1:end/2) + positionsMs{i_filename}(indLogicals{i_filename});
    node = transpose(reshape(positions,2,[]));
    
    TR = triangulation(elems{i_filename}, node(:,1), node(:,2));
    boundaries{i_filename} = freeBoundary(TR)';
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
            if strcmp(simulation_type_list{i_filename},'CG')
                u = trajectories{i_filename}(:,ti);
                
                indLogical = indLogicals{i_filename};
                positions = positionss{i_filename};
                positionsM = positionsMs{i_filename};
                positions(indLogical) = u(1:end/2) + positionsM(indLogical);
                node = transpose(reshape(positions,2,[]));
                
                trimesh(elems{i_filename},node(:,2),node(:,1),i_filename*ones(size(node,1),1),'FaceAlpha',0.0,'EdgeAlpha',0.5,'EdgeColor',edge_colors{i_filename});
                
            else % draw for DG
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
                
            end
            
            axis equal
            axis(axis_box)
            % draw a horizontal reference line for droop test
            %                 hl = refline(0,lowests{i_filename});
            %                 hl.Color = edge_colors{i_filename};
            
            drawnow
            hold on
            
        end
        legend(ha,simulation_type_list)
        frame = getframe(gca);
        writeVideo(vid,frame);
        
        
    end
    
    hold off
    
end

end