function triangle_simulation_sub

clear all
close all

draw = true;
rerun_flag = true;
save_state = true;
draw = true;

dt = 1/120;
tsteps = 120*3;

fs = filesep;
mass_flag = true;

type = 'triangle';

axis_box = [-1 1.5 -0.5 1];
% options = optimoptions('fsolve','Algorithm','levenberg-marquardt');
options = optimoptions('fsolve','TolFun',1.e-9,'TolX',1.e-9,'Display','final');

maxA = 0.1;
%         meshname = sprintf('simData%ctrimeshes%c%s maxA%.e',fs,fs,type, maxA);

% meshname = sprintf('simData%ctrimeshes%c%s maxA%.e tritool%s',fs,fs,type, maxA);

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
filename = [filename '_linear_sub'];
if (exist([filename '.mat'], 'file') ~= 2) || rerun_flag
    
    % construct triangular mesh object
    obj = elasticTriObj(nodeM, elem);
    
    
    obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neohookean
    
    obj.initGhostPoints;
%     ha = obj.init_vis;
%     obj.simple_vis(ha);
%     
%     obj.simple_vis_sub(ha);
%     obj.subElasticForce;
    M = obj.M;
    % print out the total mass for sanity check
    % should be the same for all mesh resolution
    % disp('total mass')
    % disp(sum(spdiags(M)))
    
    
    
    
    %         M = M(~ind_fix,~ind_fix); % extract the non-fixed part
    %
    Dx = 0*rand(2*N,1); % displacement field. set zero for rest state
    obj.SetCurrentState(Dx);
    %
    K = obj.StiffnessMatrix;
    %         K = K(~ind_fix,~ind_fix); % extract the non-fixed part
    
    
    
    [V,D] = eig(full(K),full(M));
    %     [~,d2] = eig(Ke,Me);
    [low_eig, permutation_indices] = sort(diag(D));
    %     low_eig2 = sort(diag(d2));
    %     eig_list(ind_N) = low_eig(1);
    V = V(:,permutation_indices);
    firstMode = V(:,4)/2;
    
    save([filename '.mat'], 'obj','firstMode');
else
    load([filename '.mat'], 'obj');
    load([filename '.mat'], 'firstMode');
end

figure
node = obj.node;
trimesh(elem,node(:,1),node(:,2),zeros(size(node,1),1),zeros(size(node,1),1));
trimesh(elem,node(:,1)+firstMode(1:2:end),node(:,2)+firstMode(2:2:end),zeros(size(node,1),1),zeros(size(node,1),1));
axis(axis_box)
axis equal
drawnow

positionsM = node';
positionsM = positionsM(:);
positions = positionsM;
externalGravity  = zeros(size(positions));
nFixed = 0;
indLogical = true(size(positions));



Dx = firstMode;
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
    
%     % k1
%         positions(indLogical) = positionsM(indLogical) + u(1:end/2);
%         v(indLogical) = u(end/2 + 1 : end);
%         
%         obj.SetCurrentState(positions - positionsM);
% 
%         Eforce = obj.subElasticForce;
%         Eforce = Eforce(indLogical);
%         Mass = obj.M;
%         Mass = Mass(indLogical,indLogical);
% %         fExternal = Mass * externalGravity;
%         
%         f1 = Eforce;
%         
%         k1 = dt * [v(indLogical); Mass\f1];
%         
%         % k2
%         positions(indLogical) = positionsM(indLogical) + u(1:end/2) + k1(1:end/2)/2;
%         v(indLogical) = u(end/2 + 1 : end) + k1(end/2+1 : end)/2;
%         
%         obj.SetCurrentState(positions - positionsM);
%         Eforce = obj.subElasticForce;
%         Eforce = Eforce(indLogical);
%         
%         f2 = Eforce;
%         
%         k2 = dt * [v(indLogical); Mass\f2];
%         
%         % k3
%         positions(indLogical) = positionsM(indLogical) + u(1:end/2) + k2(1:end/2)/2;
%         v(indLogical) = u(end/2 + 1 : end) + k2(end/2+1 : end)/2;
%         
%         obj.SetCurrentState(positions - positionsM);
%         Eforce = obj.subElasticForce;
%         Eforce = Eforce(indLogical);
%         
%         f3 = Eforce;
%         
%         k3 = dt * [v(indLogical); Mass\f3];
%         
%         % k4
%         positions(indLogical) = positionsM(indLogical) + u(1:end/2) + k3(1:end/2);
%         v(indLogical) = u(end/2 + 1 : end) + k3(end/2+1 : end);
%         
%         obj.SetCurrentState(positions - positionsM);
%         Eforce = obj.subElasticForce;
%         Eforce = Eforce(indLogical);
%         
%         f4 = Eforce;
%         
%         
%         k4 = dt * [v(indLogical); Mass\f4];
%         
%         % RK4
%         %         positions(indLogical) = positionsM(indLogical) + u(1:end/2) + k1(1:end/2)/6 + k2(1:end/2)/3 + k3(1:end/2)/3 + k4(1:end/2)/6;
%         %         v(indLogical) = u(end/2 + 1 : end) + k1(end/2 + 1 : end)/6 + k2(end/2 + 1 : end)/3 + k3(end/2 + 1 : end)/3 + k4(end/2 + 1 : end)/6;
%         u(1:end/2) = u(1:end/2) + k1(1:end/2)/6 + k2(1:end/2)/3 + k3(1:end/2)/3 + k4(1:end/2)/6;
%         u(end/2 + 1:end)= u(end/2 + 1 : end) + k1(end/2 + 1 : end)/6 + k2(end/2 + 1 : end)/3 + k3(end/2 + 1 : end)/3 + k4(end/2 + 1 : end)/6;
%         
% %     u = IM(dt, u, obj);
    
    
    if(draw)
        if or(mod(ti, draw_rate) == 1, draw_rate == 1)
        obj.SetCurrentState(Dx)
            node = u(1:end/2) + positionsM(indLogical);
            node = transpose(reshape(node,2,[]));
            trimesh(elem,node(:,1),node(:,2),zeros(size(node,1),1),zeros(size(node,1),1));
            %     trimesh(elem,node(:,1)+firstMode(1:2:end),node(:,2)+firstMode(2:2:end),zeros(size(node,1),1),zeros(size(node,1),1));
            axis(axis_box*3)
            axis equal
            drawnow
            hold on
            ef = obj.subElasticForce;

            quiver(node(:,1),node(:,2),ef(1:2:end),ef(2:2:end),0)
            hold off
            if save_state
                frame = getframe;
                writeVideo(vid,frame);
            end
        end
        
    end
end


fname = [filename 'Trajectory'];
save(fname)


    function out = IM(dt, u, obj)
        it = 0;
        MaxIT = 10;
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
        
        Eforce_mid = obj.subElasticForce;
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
        
        while (Dv'*Dv > 1e-6) && (it ~= MaxIT)
            
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
            
            Eforce_mid = obj.subElasticForce;
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