function rect_simulation

clear all
close all

draw = true;
rerun_flag = true;
save_state = true;
draw = true;

dt = 1/120;
tsteps = 120*3;

fs = filesep;

mesh_shape = 'rect';

constraints = 1; % types of constraint
% 1: free

deformation_mode_number = 1;

Y = 1000; % Young's modululs
P = 0.48; % Poisson ratio
rho = 1; % density
a = 0.0; % rayleigh damping
b = 0.0;

axis_box = [-0.5 .5 -3 1];

maxA = 0.1;

meshname = sprintf('mesh_data%c%s_maxA_%.d',fs,mesh_shape, maxA);

if exist([meshname '.mat'], 'file') ~= 2
    disp('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end

% options for fsolve in implicit solver
options = optimoptions('fsolve','TolFun',1.e-9,'TolX',1.e-9,'Display','final');
% options = optimoptions('fsolve','Algorithm','levenberg-marquardt');


% fix the orientation
elem(:,[1 3]) = elem(:,[3 1]);

N = size(nodeM,1);

filename = sprintf('sim_data%csimulation_%s_maxA_%.d', fs, mesh_shape,maxA);
if (exist([filename '.mat'], 'file') ~= 2) || rerun_flag
    
    % construct triangular mesh object
    obj = elasticTriObj(nodeM, elem);
    
    
    obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to linear
    %
    Dx = 0*rand(2*N,1); % displacement field. set zero for rest state
    obj.SetCurrentState(Dx);
    %
    M = obj.M;
    K = obj.StiffnessMatrix;
    %         K = K(~ind_fix,~ind_fix); % extract the non-fixed part
%     
%     [V,D] = eig(full(K),full(M));
%     [low_eig, permutation_indices] = sort(diag(D));
%     V = -V(:,permutation_indices);
%     firstMode = V(:,4)/2;
    
%     save([filename '.mat'], 'obj','V','D');
else
    load([filename '.mat']);
%     load([filename '.mat'], );
end

% deformation_scale_factor = 2;
% deformation_mode = V(:,3 + deformation_mode_number)/deformation_scale_factor;

deformation_mode = zeros(numel(obj.node),1);

node = obj.node;
positionsM = node';
positionsM = positionsM(:);
positions = positionsM;
externalGravity  = zeros(size(positions));
externalGravity(1:2:end) = -9.8;
nFixed = 0;
indLogical = true(size(positions));


Xind_top = (abs(nodeM(:,1)-max(nodeM(:,1))) < 1e-6);
nFixed = sum(Xind_top);
ind_fix = reshape(transpose(repmat(Xind_top,1,2)),[],1); % logical index for total position vector

indLogical(ind_fix) = false;

% files = dir(sprintf('sim_data%csimulation*.mat',fs));
% for i_file = 1:length(files)
%     loaded_struct_array = load(files(i_file).name);
%     fields = fieldnames(loaded_struct_array);
%     for i_fields = 1:numel(fields)
%         % don't check the iterators (begin with 'i')
%         if fields{i_fields}(1) ~= 'i'
%             eval(fields{i_fields} ==
% end


Dx = deformation_mode(indLogical);
v = zeros(length(Dx),1);

u = [Dx; v];


if save_state && draw
    
    vidname = strcat(filename,'.avi');
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
    u = IM(dt, u, obj);
    
    
    if(draw)
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)

            positions(indLogical) = u(1:end/2) + positionsM(indLogical);
            node = transpose(reshape(positions,2,[]));
            triplot(elem,node(:,2),node(:,1),zeros(size(node,1),1),zeros(size(node,1),1));
            %     trimesh(elem,node(:,1)+firstMode(1:2:end),node(:,2)+firstMode(2:2:end),zeros(size(node,1),1),zeros(size(node,1),1));
            axis(axis_box)
            axis equal
            drawnow
            if save_state
                frame = getframe(gcf);
                writeVideo(vid,frame);
            end
        end
        
    end
end


% fname = [filename '_trajectory.mat'];
% save(fname)

    save([filename '.mat']);


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
        
        obj.SetCurrentState(positions - positionsM);
        
        K_mid = obj.StiffnessMatrix;
        K_mid = 1/2 * (K_mid + K_mid');
        %         K0 = obj.restStiffness;
        Mass = obj.M;
        Mass = Mass(indLogical,indLogical);
        K_mid = K_mid(indLogical,indLogical);
        %         K0 = obj.K0;
        %         B = -a * Mass - b * K0(indLogical,indLogical);
        B = -a * Mass - b * K_mid;
        
        Eforce_mid = obj.ElasticForce;
        Eforce_mid = Eforce_mid(indLogical);
        
        fExternal = Mass * externalGravity(indLogical);
        
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
            %             v_free_mid = 1/2 * (v_free_new + v_free);
            %             v(indLogical) = v_free_mid;
            positions(indLogical) = positionsM(indLogical) + dq_free + 1/4*dt*(v_free + v_free_new);
            
            obj.SetCurrentState(positions - positionsM);
            
            K_mid = obj.StiffnessMatrix;
            K_mid = 1/2 * (K_mid + K_mid');
            Mass = obj.M;
            Mass = Mass(indLogical,indLogical);
            K_mid = K_mid(indLogical,indLogical);
            
            K0 = obj.K0;
            %             B = -a * Mass - b * K0(indLogical,indLogical);
            B = -a * Mass - b * K_mid;
            
            Eforce_mid = obj.ElasticForce;
            Eforce_mid = Eforce_mid(indLogical);
            
            fExternal = Mass * externalGravity(indLogical);
            
            f_mid = Eforce_mid + fExternal + B*1/2*(v_free+v_free_new);
            
            Dv = -(speye(2*(N-nFixed)) + 1/4* dt*dt*(Mass\K_mid) - 1/2 * dt*(Mass\B))\(v_free_new - v_free - dt * (Mass\f_mid));
            v_free_new = v_free_new + Dv;
            
            residual = (v_free_new - v_free - dt * (Mass\f_mid))' * (v_free_new - v_free - dt * (Mass\f_mid));
            %             residual_list = [residual_list residual];
            %             Dv_norm_list = [Dv_norm_list Dv'*Dv];
            it = it + 1
            
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