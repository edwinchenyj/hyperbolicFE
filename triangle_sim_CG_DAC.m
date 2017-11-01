function triangle_sim_CG_DAC
%% set parameters, flags, and meshes
close all
fs = filesep;

%% load the eigenvalues from the fine mesh
maxA = 0.001; % resolution of the fine mesh
mesh_shape = 'triangle';
ela_material = 'linear'; % choice: 'linear', 'neo-hookean'
dirname = sprintf('sim_data%c%s_%s_maxA_%.d', fs, ela_material, mesh_shape, maxA);
load([dirname fs 'data.mat']);
low_eig_fine = low_eig;

%%

% draw = false;
rerun_flag = true;
save_state = true;
draw = true;
% test_mode = true;

dt = 1/120;
tsteps = 120*3;

fs = filesep;

mesh_shape = 'triangle';
maxA = 0.001;
simulation_type = 'CG';

% number of eigenvalue fitted
NE = 1;
DAC_on = true;

solver = 'SI';
constraints = 1; % types of constraint
% 1: free


Y = 100; % Young's modululs
P = 0.48; % Poisson ratio
rho = 1; % density
a = 0.0; % rayleigh damping
b = 0.00;
ela_material = 'linear'; % choice: 'linear', 'neo-hookean'

axis_box = [-1 1.5 -0.5 1];

meshname = sprintf('mesh_data%c%s_maxA_%.d',fs,mesh_shape, maxA);

if exist([meshname '.mat'], 'file') ~= 2
    disp('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end


%%

N = size(nodeM,1);

dirname = sprintf('sim_data%c%s_%s_maxA_%.d', fs, ela_material, mesh_shape, maxA);
if (exist([dirname fs 'data.mat'], 'file') ~= 2) || rerun_flag
    
    % construct triangular mesh object
    obj = elasticTriObj(nodeM, elem);
    switch ela_material
        case 'linear'
            obj.SetMaterial( Y, P, rho, 1:size(elem,1), 2); % set the tri to linear
            %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
        case 'neo-hookean'
            obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
    end
    %
    Dx = 0*rand(2*N,1); % displacement field. set zero for rest state
    obj.SetCurrentState(Dx);
    %
    M = obj.M;
    K = obj.StiffnessMatrix;
    
    [V,D] = eig(full(K),full(M));
    [low_eig, permutation_indices] = sort(diag(D));
    V = -V(:,permutation_indices);
    
    
    mkdir(dirname);
    save([dirname fs 'data.mat'], 'obj','V','low_eig'); % storing eigen decomp
else
    load([dirname fs 'data.mat']);
end

if maxA ~= 0.001
p = polyfit(low_eig(4:NE+3),low_eig_fine(4:NE+3),NE-1);
end
ha = obj.init_vis;



% deformation for the initial condition
deformation_mode = zeros(size(obj.X));
for di = 1:3
deformation_mode = deformation_mode + V(:,3 + di)/di; % divide the mode by it's mode number to reduce the amplitude of the high modes
end

if maxA == 0.1
    deformation_mode = -deformation_mode;
end

Dx = deformation_mode;


node = obj.node;

positionsM = node';
positionsM = positionsM(:);
positions = positionsM;
externalGravity  = zeros(size(positions));

indLogical = true(size(positions)); % TODO: need to map the indLogical to the indLogical for DG
nFixed = 0;

obj.current_vis(obj.vis_handle);

v = zeros(length(Dx),1);
axis(axis_box)
axis equal
xlim = ha.XLim;
ylim = ha.YLim;
u = [Dx; v];

frame_rate = 60;
if save_state && draw
    if DAC_on
    simdir = strcat(dirname,fs,solver,'_',simulation_type,'_DAC',num2str(NE),'_Y',num2str(Y),'_P',num2str(P),'_dt',num2str(dt));
    else
        simdir = strcat(dirname,fs,solver,'_',simulation_type,'_DAC_off_Y',num2str(Y),'_P',num2str(P),'_dt',num2str(dt));
    end
    mkdir(simdir);
    vidname = strcat(simdir,fs,'video.avi');
    vid = VideoWriter(vidname);
    
    vid.FrameRate = frame_rate;
    open(vid);
end

% rate to draw the scene
sim_rate = round(1/dt);
draw_rate = round(sim_rate/frame_rate);

trajectory = zeros(size(u,1),tsteps);
for ti = 1:tsteps
    trajectory(:,ti) = u;
    
    switch solver
        case 'IM'
            u = IM(dt, u, obj);
        case 'SI'
            if strcmp(ela_material,'neo-hookean')
                error('neo-hookean not supported now');
            else
                if (ti == 1 && DAC_on)
                    if maxA ~= 0.001
                    K = obj.StiffnessMatrix;
                    Minv_K_new = polyvalm(p,M\K);
                    [v_new, d_new] = eig(full(Minv_K_new),full(M));
                    K = M * Minv_K_new;
                    else
                    K = obj.StiffnessMatrix;
                    end
                end
                Eforce = -K*(obj.x-obj.X);
                Mass = M(indLogical,indLogical);
                K = K(indLogical,indLogical);
                B = -a * Mass - b * K;
                
                Eforce = Eforce(indLogical);
                
                fExternal = Mass * externalGravity;
                
                
                f = Eforce + fExternal + B*u(end/2 + 1:end); % from column to row
                
                A = (Mass - dt * B + dt^2 * K);
                v_free = u(end/2 + 1:end);
                rhs = dt * (f - dt * K * v_free);
                dv_free = A\rhs;
                
                u(end/2 +1 :end) = u(end/2 +1 :end) + dv_free;
                u(1:end/2) = u(1:end/2) + dt * u(end/2+1:end);
                dq_free = u(1:end/2);
                v_free = u(end/2+1:end);
                u_new = u;
                dq_free_new = u_new(1:end/2);
                v_free_new = u_new(end/2+1:end);
                
                positions(indLogical) = positionsM(indLogical) + dq_free + 1/4 * dt * (v_free + v_free_new);
                
                obj.SetCurrentState(positions - positionsM);
                
            end
        case 'SIIMEX'
            K = obj.StiffnessMatrix;
            Mass = obj.M;
            Eforce = obj.ElasticForce;
            
            Mass = Mass(indLogical,indLogical);
            K = K(indLogical,indLogical);
            B = -a * Mass - b * K;
            
            Eforce = Eforce(indLogical);
            
            fExternal = Mass * externalGravity;
            
            
            f = Eforce + fExternal + B*u(end/2 + 1:end); % from column to row
            
            A = (Mass - dt * B + dt^2 * K);
            v_free = u(end/2 + 1:end);
            rhs = dt * (f - dt * K * v_free);
            dv_free = A\rhs;
            
            u(end/2 +1 :end) = u(end/2 +1 :end) + dv_free;
            u(1:end/2) = u(1:end/2) + dt * u(end/2+1:end);
            dq_free = u(1:end/2);
            v_free = u(end/2+1:end);
            u_new = u;
            dq_free_new = u_new(1:end/2);
            v_free_new = u_new(end/2+1:end);
            
            positions(indLogical) = positionsM(indLogical) + dq_free + 1/4 * dt * (v_free + v_free_new);
            
            obj.SetCurrentState(positions - positionsM);
    end
    
    if(draw)
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)
            axis(axis_box)
            cla
            obj.current_vis(obj.vis_handle);
            obj.vis_handle.XLim = xlim;
            obj.vis_handle.YLim = ylim;
            if save_state
                frame = getframe(gca);
                writeVideo(vid,frame);
            else
                drawnow;
            end
        end
        
    end
end


if save_state
    save([simdir fs 'trajectory.mat']);
end



end