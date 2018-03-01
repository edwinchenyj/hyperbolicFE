function filename_list=bar_error

clear all
close all

pname = 'bar_error';

fs = filesep;

axis_box = [-1 1 -1 1 -1 1]/4;
% compare 2 simulations
filename_list = cell(1,6);

%% first comparision
filename_list{1} = 'sim_data/simulation3D_deform_full_SI_not_rescaled_bar26_neo-hookean_small_bar/SI_CG_constraint_n_h0.5_Y10000_P0.45_rho1000_a0_b0_dt0.01_def-scln_redu20/trajectory.mat';

%% first comparision

filename_list{2} =     'sim_data/simulation3D_deform_full_SI_not_rescaled_bar26_neo-hookean_small_bar/SI_CG_constraint_n_h0.3_Y10000_P0.45_rho1000_a0_b0_dt0.01_def-scln_redu20/trajectory.mat';


%% second comparison

filename_list{3} = 'sim_data/simulation3D_deform_full_SI_not_rescaled_bar26_neo-hookean_small_bar/SI_CG_constraint_n_h0.2_Y10000_P0.45_rho1000_a0_b0_dt0.01_def-scln_redu20/trajectory.mat';


%% first comparision
filename_list{4} = 'sim_data/simulation3D_deform_full_SI_rescaled_bar220_neo-hookean_small_bar/SI_CG_constraint_n_h0.48_Y50000_P0.45_rho1000_a0_b0_dt0.01_def-scln_redu20/trajectory.mat';

%% first comparision

filename_list{5} = 'sim_data/simulation3D_deform_full_SI_rescaled_bar220_neo-hookean_small_bar/SI_CG_constraint_n_h0.3_Y50000_P0.45_rho1000_a0_b0_dt0.01_def-scln_redu20/trajectory.mat';


%% second comparison

filename_list{6} = 'sim_data/simulation3D_deform_full_SI_rescaled_bar220_neo-hookean_small_bar/SI_CG_constraint_n_h0.2_Y50000_P0.45_rho1000_a0_b0_dt0.01_def-scln_redu20/trajectory.mat';



%%
colors = [summer(3); bone(3); copper(3)];
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
        disp('NT')
        S.obj.NT
        disp('N')
        S.obj.N
    end
    close(gcf)
    min_positions_ind = zeros(size(filename_list));
    % boundary for the mesh in CG sim
    boundaries = cell(size(filename_list));
    for i_filename = 1:length(trajectories)
        u = trajectories{i_filename}(:,100);
        indLogical = indLogicals{i_filename};
        positions = positionss{i_filename};
        
        positions = u(1:end/2) + positionsMs{i_filename};
        node = transpose(reshape(positions,3,[]));
        [~,ind] = min(positions(1:3:end));
        min_positions_ind(i_filename) = ind;
%         TR = triangulation(elems{i_filename}, node(:,1), node(:,2));
%         boundaries{i_filename} = freeBoundary(TR)';
    end
    
    
    tsteps = 200;
    minX_traj = zeros(tsteps,length(trajectories));
    for ti = 1:tsteps
        
        
            for i_filename = 1:length(trajectories)

                u = trajectories{i_filename}(:,ti);
                
                indLogical = indLogicals{i_filename};
                positions = positionss{i_filename};
                positionsM = positionsMs{i_filename};
                positions = u(1:end/2) + positionsM;
                minX_traj(ti,i_filename) = positions(3*(min_positions_ind(i_filename)-1)+1);
            end
                
    end
    xscale = max(positionsM(1:3:end)) - min(positionsM(1:3:end));
    
    close all
    minX_traj(:,1:2) = minX_traj(:,1:2) - repmat(minX_traj(:,3),1,2);
    minX_traj(:,4:5) = minX_traj(:,4:5) - repmat(minX_traj(:,6),1,2);
    minX_traj = minX_traj(:,[1 2 4 5]);
    plot(abs(minX_traj/xscale),'*-','linewidth',5)
    disp('max err')
    max(abs(minX_traj/xscale))
    ah = gca;
    ah.FontSize = 20;
    xlabel('Frame #')
    ylabel('Error');
%     max(maxZ_traj)
    legend({'plain FEM 159 ','plain FEM 572','NMSC 159','NMSC 572'},'FontSize',24)
        print(pname,'-dpdf')        
        print(pname,'-dpng')
        savefig(pname)
    
    
end
    
