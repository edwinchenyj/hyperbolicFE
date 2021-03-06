function four_elem_rev_sim
%% set parameters, flags, and meshes
close all

rerun_flag = true;
save_state = true;
draw = true;

dt = 1/120;
tsteps = 120*1;

fs = filesep;

mesh_shape = 'four_elem_rev';
simulation_type = 'CG';

% set the DG flag base on simulation type
switch simulation_type(1:2)
    case 'DG'
        isDG = true;
    otherwise
        isDG = false;
end

DGeta = 1e2;
solver = 'SI';
constraints = 1; % types of constraint
% 1: free
deformation_scale_factor = 10;
deformation_mode_number = 1;

Y = 100; % Young's modululs
P = 0.48; % Poisson ratio
rho = 1; % density
a = 0.0; % rayleigh damping
b = 0.00;
material = 'linear'; % choice: 'linear', 'neo-hookean'

axis_box = [-1 1.5 -0.5 1];

elem = [1 2 3; 2 4 3; 2 6 4; 2 5 6];
nodeM = [0 0; 1 0; 0 1; 1 1; 2 0; 2 1]/2 ;
% nodeM(:,1) = nodeM(:,1) - 1/2;
% nodeM(:,2) = nodeM(:,2) - 1/4;
%%


% fix the orientation
% elem(:,[1 3]) = elem(:,[3 1]);

N = size(nodeM,1);

dirname = sprintf('sim_data%c%s_%s', fs, material, mesh_shape);
if (exist([dirname fs 'data.mat'], 'file') ~= 2) || rerun_flag
    
    % construct triangular mesh object
    if isDG
        obj = elasticDGTriObj(nodeM, elem);
        obj.eta = DGeta;
    else
        obj = elasticTriObj(nodeM, elem);
    end
    switch material
        case 'linear'
            obj.SetMaterial( Y, P, rho, 2, 0, 0); % set the tri to linear
            %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
        case 'neo-hookean'
            obj.SetMaterial( Y, P, rho, 1, 0, 0); % set the tri to neo-hookean
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
    firstMode = V(:,4)/2;
    
    mkdir(dirname);
    save([dirname fs 'data.mat'], 'obj','V','D'); % storing eigen decomp
else
    load([dirname fs 'data.mat']);
    %     load([filename '.mat'], );
end

ha = obj.init_vis;



% deformation for the initial condition

deformation_mode = V(:,3 + deformation_mode_number)/deformation_scale_factor;

% positionsM = node';
% positionsM = positionsM(:);
% positions = positionsM;
% externalGravity  = zeros(size(positions));
% nFixed = 0;
% indLogical = true(size(positions));


Dx = deformation_mode;

if isDG
    N = size(obj.DGnodeM,1);
    node = obj.DGnode;
else
    
    node = obj.node;
end

positionsM = node';
positionsM = positionsM(:);
positions = positionsM;
externalGravity  = zeros(size(positions));

if isDG
    indLogical = true(size(positions)); % TODO: need to map the indLogical to the indLogical for DG
    nFixed = 0;
    
    Dx = obj.CGxToDGx(Dx);
    obj.SetCurrentDGState(Dx);
    obj.DG_current_vis(obj.vis_handle);
else
    indLogical = true(size(positions)); % TODO: need to map the indLogical to the indLogical for DG
    nFixed = 0;
    obj.SetCurrentState(Dx);
    obj.simple_vis(obj.vis_handle);
    
end

v = zeros(length(Dx),1);
axis(axis_box)
axis equal
xlim = ha.XLim;
ylim = ha.YLim;
u = [Dx; v];

if save_state && draw
    if isDG
        simdir = strcat(dirname,fs,solver,'_',simulation_type,'_Y',num2str(Y),'_P',num2str(P),'_dt',num2str(dt),'_eta',num2str(DGeta));
    else
        simdir = strcat(dirname,fs,solver,'_',simulation_type,'_Y',num2str(Y),'_P',num2str(P),'_dt',num2str(dt));
    end
    mkdir(simdir);
    vidname = strcat(simdir,fs,'video.avi');
    vid = VideoWriter(vidname);
    vid.FrameRate = 60;
    open(vid);
end

% rate to draw the scene
sim_rate = round(1/dt);
draw_rate = round(sim_rate/vid.FrameRate);

trajectory = zeros(size(u,1),tsteps);
for ti = 1:tsteps
    trajectory(:,ti) = u;
    
    switch solver
        case 'IM'
            u = IM(dt, u, obj);
        case 'SI'
            if isDG
                K = obj.DGStiffnessMatrix;
                Mass = obj.DGM;
                Eforce = obj.DGElasticForce;
            else
                K = obj.StiffnessMatrix;
                Mass = obj.M;
                Eforce = obj.ElasticForce;
            end
            
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
            
            if isDG
                obj.SetCurrentDGState(positions - positionsM);
            else
                obj.SetCurrentState(positions - positionsM);
                
            end
            
    end
    if(draw)
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)
            axis(axis_box)
            cla
            if isDG
                obj.DG_current_vis(obj.vis_handle);
            else
                obj.simple_vis(obj.vis_handle);
            end
            obj.vis_handle.XLim = xlim;
            obj.vis_handle.YLim = ylim;
            if save_state
                frame = getframe(gca);
                writeVideo(vid,frame);
            end
        end
        
    end
end


% fname = [filename '_trajectory.mat'];
% save(fname)

save([simdir fs 'trajectory.mat']);


% integrators
    function out = IM(dt, u, obj)
        it = 0;
        MaxIT = 40;
        Dv = Inf;
        dq_free = u(1:end/2);
        v_free = u(end/2+1:end);
        u_new = u;
        dq_free_new = u_new(1:end/2);
        v_free_new = u_new(end/2+1:end);
        
        positions(indLogical) = positionsM(indLogical) + dq_free + 1/4 * dt * (v_free + v_free_new);
        
        
        if isDG
            obj.SetCurrentDGState(positions - positionsM);
            K_mid = obj.DGStiffnessMatrix;
            Mass = obj.DGM;
            Eforce_mid = obj.DGElasticForce;
        else
            obj.SetCurrentState(positions - positionsM);
            K_mid = obj.StiffnessMatrix;
            Mass = obj.M;
            Eforce_mid = obj.ElasticForce;
        end
        
        Mass = Mass(indLogical,indLogical);
        K_mid = K_mid(indLogical,indLogical);
        B = -a * Mass - b * K_mid;
        
        Eforce_mid = Eforce_mid(indLogical);
        
        fExternal = Mass * externalGravity;
        
        f_mid = Eforce_mid + fExternal + B*1/2*(v_free + v_free_new);
        
        residual0 = (dt * (Mass\f_mid))' * (dt * (Mass\f_mid));
        Dv = -(speye(2*(N-nFixed)) + 1/4* dt*dt*(Mass\K_mid) - 1/2 * dt*(Mass\B))\(v_free_new - v_free - dt * (Mass\f_mid));
        v_free_new = v_free_new + Dv;
        
        residual = (v_free_new - v_free - dt * (Mass\f_mid))' * (v_free_new - v_free - dt * (Mass\f_mid));
        
        it = it + 1;
        
        u_new(1:end/2) = dq_free + 1/2 * dt * (v_free + v_free_new);
        u_new(end/2+1:end) = v_free_new;
        
        while (Dv'*Dv > 1e-12) && (it ~= MaxIT)
            
            v_free_new = u_new(end/2+1:end);
            positions(indLogical) = positionsM(indLogical) + dq_free + 1/4*dt*(v_free + v_free_new);
            
            if isDG
                obj.SetCurrentDGState(positions - positionsM);
            else
                obj.SetCurrentState(positions - positionsM);
            end
            
            if isDG
                K_mid = obj.DGStiffnessMatrix;
                Mass = obj.DGM;
                Eforce_mid = obj.DGElasticForce;
            else
                K_mid = obj.StiffnessMatrix;
                Mass = obj.M;
                Eforce_mid = obj.ElasticForce;
            end
            
            Mass = Mass(indLogical,indLogical);
            K_mid = K_mid(indLogical,indLogical);
            
            B = -a * Mass - b * K_mid;
            
            Eforce_mid = Eforce_mid(indLogical);
            
            fExternal = Mass * externalGravity;
            
            f_mid = Eforce_mid + fExternal + B*1/2*(v_free+v_free_new);
            
            Dv = -(speye(2*(N-nFixed)) + 1/4* dt*dt*(Mass\K_mid) - 1/2 * dt*(Mass\B))\(v_free_new - v_free - dt * (Mass\f_mid));
            v_free_new = v_free_new + Dv;
            
            residual = (v_free_new - v_free - dt * (Mass\f_mid))' * (v_free_new - v_free - dt * (Mass\f_mid));
            %             residual_list = [residual_list residual];
            %             Dv_norm_list = [Dv_norm_list Dv'*Dv];
            it = it + 1;
            
            u_new(1:end/2) = dq_free + 1/2 * dt * (v_free + v_free_new);
            u_new(end/2+1:end) = v_free_new;
            
            if (it > 3 && residual > residual0) || it == MaxIT
                disp('local substep required')
                u_half = IM(dt/2, u, obj);
                
                v_free = u_half(end/2 + 1:end);
                dq_free = u_half(1:end/2);
                
                v(indLogical) = v_free;
                positions(indLogical) = positionsM(indLogical) + dq_free;
                
                u_half = [positions(indLogical)-positionsM(indLogical); v(indLogical)];
                
                u_new = IM(dt/2, u_half, obj);
                break;
            end
            
            if it == MaxIT
                disp('Newton iteration not converging in IM')
                
                break;
            end
            
        end
        
        out = u_new;
    end


end