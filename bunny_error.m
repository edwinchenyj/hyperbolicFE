function filename_list=bunny_error

clear all
close all

pname = 'bunny_error';

fs = filesep;

axis_box = [-1 1 -1 1 -1 1]/4;
% compare 2 simulations
filename_list = cell(1,4);


%% first comparision

filename_list{1} = 'sim_data/simulation3D_jump_not_reduced_SI_vcoarsebunny_stvk/SI_CG_constraint_n_h0.2_Y700_P0.45_rho1000_a0_b0_dt0.01_def-scln/trajectory.mat';



%% second comparison

filename_list{2} = 'sim_data/simulation3D_jump_SI_rescaled_vcoarsebunny_stvk20/SI_CG_constraint_n_h0.2_Y700_P0.45_rho1000_a0_b0_dt0.01_def-scln/trajectory.mat';

%% second comparison

filename_list{3} = 'sim_data/simulation3D_jump_SI_DAC_vcoarsebunny_stvk20/SI_CG_constraint_n_h0.2_Y2024.976_P0.45_rho1000_a0_b0_dt0.01_def-scln/trajectory.mat';

%% second comparison

filename_list{4} = 'sim_data/simulation3D_jump_not_reduced_SI_bunny_stvk/SI_CG_constraint_n_h0.2_Y700_P0.45_rho1000_a0_b0_dt0.01_def-scln/trajectory.mat';

%%
colors = [summer(3); autumn(3); winter(3)];
colors = colors(2:3:end,:);
edge_colors = num2cell(colors,2); % create a list of different colours
trajectories = cell(size(filename_list));
positionsMs = cell(size(filename_list));
indLogicals = cell(size(filename_list));
elems = cell(size(filename_list)); % faces
positionss = cell(size(filename_list));
axis_boxes = cell(size(filename_list));
lowests= cell(size(filename_list));
h = cell(size(filename_list));



 
    for i_filename = 1:length(filename_list)
        filename = filename_list{i_filename};
        
        S = load(filename);
        cla
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
    close(gcf)
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
    
    
    tsteps = 30;
    maxZ_traj = zeros(tsteps,length(trajectories));
    for ti = 1:tsteps
        
        
            for i_filename = 1:length(trajectories)

                u = trajectories{i_filename}(:,ti);
                
                indLogical = indLogicals{i_filename};
                positions = positionss{i_filename};
                positionsM = positionsMs{i_filename};
                positions = u(1:end/2) + positionsM;
                maxZ_traj(ti,i_filename) = max(positions(3:3:end));
            end
                
    end
        zscale = max(positionsM(3:3:end)) - min(positionsM(3:3:end));

    close all
    maxZ_traj = maxZ_traj - repmat(maxZ_traj(:,end),1,4);
    maxZ_traj = maxZ_traj(:,1:2);
    plot(abs(maxZ_traj/zscale),'*-','linewidth',5)
    disp('max err')
    max(abs(maxZ_traj/zscale))
    
    ah = gca;
    ah.FontSize = 20;
    xlabel('Frame #')
    ylabel('Error');
%     max(maxZ_traj)
    legend({'plain FEM','NMSC'},'FontSize',24)
        print(pname,'-dpdf')        
        print(pname,'-dpng')
        savefig(pname)
    
    
end
    
