function filename_list=video_comparsion_script

clear all
close all

video_name = 'CG_polyfit_1_2_mode_comparison_gravity2_P48.avi';

fs = filesep;

% compare four simulations
filename_list = cell(1,4);

%% first comparision
draw = true;
rerun_flag = true;
save_state = true;
dt = [1/100];
T = [2];
mesh_shape = 'rect_horizontal';
simulation_type = 'CG';
solver = 'IM';

Y = [6e4];
P = [0.48];
rho = [1000];
a = [0.01];
b = [0.005];
constraint = 'top';
material_type = 'neo-hookean';
gravity = 'on';
DGeta = 1e6;
axis_box = [-0.5 .5 -2 1];
deformation_scale = 10;


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

maxA = [0.001];

dirname = sprintf('sim_data%c%s_%s', fs, material_type, mesh_shape);
if isDG
    simdir = strcat(dirname,fs,solver,'_',...
        simulation_type,...
        '_constraint_',constraint,...
        '_maxA',num2str(maxA),...
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
        '_maxA',num2str(maxA),...
        '_Y',num2str(Y),...
        '_P',num2str(P),...
        '_rho',num2str(rho),...
        '_a',num2str(a),...
        '_b',num2str(b),...
        '_dt',num2str(dt),...
        '_def-scl',num2str(deformation_scale));
    
end
filename_list{1} = [simdir fs 'trajectory.mat'];


%% second comparision


simulation_type = 'CG';
solver = 'IM';
maxA = [0.1];
polyfit_modes = 1;
is_polyfit = true;
dirname = sprintf('sim_data%c%s_%s', fs, material_type, mesh_shape);


if isDG
    simdir = strcat(dirname,fs,solver,'_',...
        '_polyfit_',...
        simulation_type,...
        '_constraint_',constraint,...
        '_maxA',num2str(maxA),...
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
        '_polyfit_',...
        '_constraint_',constraint,...
        '_maxA',num2str(maxA),...
        '_Y',num2str(Y),...
        '_P',num2str(P),...
        '_rho',num2str(rho),...
        '_a',num2str(a),...
        '_b',num2str(b),...
        '_dt',num2str(dt),...
        '_def-scl',num2str(deformation_scale),...
        'polyfit_modes',num2str(polyfit_modes));
    
end

filename_list{2} = [simdir fs 'trajectory.mat'];


%% third comparison

simulation_type = 'CG';
solver = 'IM';
maxA = [0.1];
polyfit_modes = 2;
is_polyfit = true;
dirname = sprintf('sim_data%c%s_%s', fs, material_type, mesh_shape);


if isDG
    simdir = strcat(dirname,fs,solver,'_',...
        '_polyfit_',...
        simulation_type,...
        '_constraint_',constraint,...
        '_maxA',num2str(maxA),...
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
        '_polyfit_',...
        '_constraint_',constraint,...
        '_maxA',num2str(maxA),...
        '_Y',num2str(Y),...
        '_P',num2str(P),...
        '_rho',num2str(rho),...
        '_a',num2str(a),...
        '_b',num2str(b),...
        '_dt',num2str(dt),...
        '_def-scl',num2str(deformation_scale),...
        'polyfit_modes',num2str(polyfit_modes));
    
end

filename_list{3} = [simdir fs 'trajectory.mat'];

%% fourth comparision

simulation_type = 'CG';
solver = 'IM';
maxA = [0.1];

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


dirname = sprintf('sim_data%c%s_%s', fs, material_type, mesh_shape);
if isDG
    simdir = strcat(dirname,fs,solver,'_',...
        simulation_type,...
        '_constraint_',constraint,...
        '_maxA',num2str(maxA),...
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
        '_maxA',num2str(maxA),...
        '_Y',num2str(Y),...
        '_P',num2str(P),...
        '_rho',num2str(rho),...
        '_a',num2str(a),...
        '_b',num2str(b),...
        '_dt',num2str(dt),...
        '_def-scl',num2str(deformation_scale));
    
end
filename_list{4} = [simdir fs 'trajectory.mat'];
%%

edge_colors = num2cell(hsv(4),2); % create a list of different colours
trajectories = cell(size(filename_list));
positionsMs = cell(size(filename_list));
indLogicals = cell(size(filename_list));
elems = cell(size(filename_list)); % faces
positionss = cell(size(filename_list));
axis_boxes = cell(size(filename_list));
lowests= cell(size(filename_list));


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
        node = transpose(reshape(positions,2,[]));
        
%         TR = triangulation(elems{i_filename}, node(:,1), node(:,2));
%         boundaries{i_filename} = freeBoundary(TR)';
    end
    
    
    % rate to draw the scene
    sim_rate = round(1/dt);
    draw_rate = round(sim_rate/vid.FrameRate);
    axis_box = axis_boxes{1};
    axis_box(4) = 1;
    axis_box(1) = -0.5;
    axis_box(2) = 1.5;
    
    hf = figure;
    
    for ti = 1:tsteps
        
        
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)
            
            for i_filename = 1:length(trajectories)
                u = trajectories{i_filename}(:,ti);
                
                indLogical = indLogicals{i_filename};
                positions = positionss{i_filename};
                positionsM = positionsMs{i_filename};
                positions = u(1:end/2) + positionsM;
                node = transpose(reshape(positions,2,[]));
                
                triplot(elems{i_filename},node(:,1),node(:,2),'Color',edge_colors{i_filename});
                
                axis equal
                axis(axis_box)
                % draw a horizontal reference line for droop test
                %                 hl = refline(0,lowests{i_filename});
                %                 hl.Color = edge_colors{i_filename};
                
                drawnow
                hold on
            end
            
            
            legend('fine','coarse with 1 polyfit','coarse with 2 polyfit','coarse','Location','southeast')
            frame = getframe(gca);
            writeVideo(vid,frame);
        end
        
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