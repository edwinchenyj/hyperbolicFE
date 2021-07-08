function triangle_sim_siere
%% set parameters, flags, and meshes
close all

% draw = false;
rerun_flag = true;
save_state = true;
draw = true;
% test_mode = true;

dt = 1/100;
tsteps = 1000*3;

fs = filesep;

mesh_shape = 'triangle';
maxA = 0.01;
simulation_type = 'CG';

solver = 'SIERE';
modes = 3;
% constraints = 1; % types of constraint
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
Y = 1e5; % Young's modululs
P = 0.45; % Poisson ratio
rho = 1000; % density
a = 0.000; % rayleigh damping
b = 0.005;
material = 'neo-hookean'; % choice: 'linear', 'neo-hookean'

axis_box = [-1 1.5 -0.5 1];

meshname = sprintf('mesh_data%c%s_maxA_%.d',fs,mesh_shape, maxA);

if exist([meshname '.mat'], 'file') ~= 2
    disp('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end


%%

N = size(nodeM,1);

dirname = sprintf('sim_data%c%s_%s', fs, material, mesh_shape);
if (exist([dirname fs 'data.mat'], 'file') ~= 2) || rerun_flag
    
    % first construct the CG object for eigen decomps
    obj = elasticTriObj(nodeM, elem);
    switch material
        case 'linear'
            obj.SetMaterial( Y, P, rho, 2, a, b); % set the tri to linear
            %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
        case 'neo-hookean'
            obj.SetMaterial( Y, P, rho, 1, a, b); % set the tri to neo-hookean
    end
    %
    Dx = 0*rand(2*N,1); % displacement field. set zero for rest state
    obj.SetCurrentState(Dx);
    %
    M = obj.M;
    K = obj.StiffnessMatrix;
    %         K = K(~ind_fix,~ind_fix); % extract the non-fixed part
    
    [V,D] = eig(full(K),full(M));
    [low_eig, permutation_indices] = sort(diag(D));
    V = -V(:,permutation_indices);
%     firstMode = V(:,4)/2;
    
    mkdir(dirname);
    save([dirname fs 'data.mat'], 'obj','V','D'); % storing eigen decomp
else
    load([dirname fs 'data.mat']);
    %     load([filename '.mat'], );
end

ha = obj.init_vis;

constraint_indices = repmat(nodeM(:,1) < 0.001,1,2);
constraint_indices = reshape((constraint_indices'),[],1);
% deformation for the initial condition

deformation_mode= zeros(size(V,1),1);

for num = deformation_mode_number
    deformation_mode = V(:,3 + num) + deformation_mode;
end
Dx = deformation_mode;

obj.SetCurrentState(Dx);
obj.simple_vis(obj.vis_handle);


v = zeros(length(Dx),1);
axis(axis_box)
axis equal
xlim = ha.XLim;
ylim = ha.YLim;
u = [Dx; v];

constraint = 'none';

if save_state && draw
    
        simdir = strcat(dirname,fs,solver,'_',simulation_type,...
            '_constraint_',constraint,...
            '_maxA',num2str(maxA),...
            '_Y',num2str(Y),...
            '_P',num2str(P),...
            '_rho',num2str(rho),...
            '_a',num2str(a),...
            '_b',num2str(b),...
            '_dt',num2str(dt));
    
    mkdir(simdir);
    vidname = strcat(simdir,fs,'video.avi');
    vid = VideoWriter(vidname);
    vid.FrameRate = 50;
    open(vid);
end


% rate to draw the scene
sim_rate = round(1/dt);


trajectory = zeros(size(u,1),tsteps);
ritz_errors = zeros(modes,tsteps);
recompute = true;
recompute_count = 1;
for ti = 1:tsteps
    trajectory(:,ti) = u;
    
    switch solver
        case 'ERE'
            u = ERE(dt, u, obj);
        case 'SIERE'
            
               u = SIERE(dt, u, obj, modes, constraint_indices,recompute);
               ritz_errors(:,ti) = obj.ritz_errors';
               if norm(obj.ritz_errors) > 1e3
                   recompute = true;
                   recompute_count = recompute_count + 1;
                   disp(recompute_count)
               else
                   recompute = false;
               end
    end
    if(draw)
        draw_rate = round(sim_rate/vid.FrameRate);
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)
            axis(axis_box)
            cla
            obj.simple_vis(obj.vis_handle);
            obj.vis_handle.XLim = xlim;
            obj.vis_handle.YLim = ylim;
            if save_state
                frame = getframe(gca);
                writeVideo(vid,frame);
            end
        end
        
    end
end
figure
plot(1:tsteps,ritz_errors');
disp(recompute_count)

% fname = [filename '_trajectory.mat'];
% save(fname)

% save([simdir fs 'trajectory.mat']);



end