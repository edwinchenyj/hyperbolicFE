function u_new = reduced_ERE_rescaled( dt, u, obj, varargin)
% inputs:
%   dt: step size
%    u: current state
%  obj: the elastic obj
%  varargin: optional input for constraints and max iteration count
% notice sum(nnz(~constraint_indices)) = length(u)/2

if nargin == 3
    constraint_indices = false(size(u,1)/2,1);
    MaxIT = 40;
elseif nargin == 4
    constraint_indices = varargin{1};
    MaxIT = 40;
elseif nargin == 5
    constraint_indices = varargin{1};
    MaxIT = varargin{2};
end

indLogical = ~constraint_indices;

Dim = obj.Dim;
reduced_dimension = obj.reduced_dimension;


dq = u(1:end/2);
v = u(end/2+1:end);

N = obj.N;
nFixed = sum(constraint_indices)/Dim;

obj.x = obj.X + dq;

obj.SetCurrentState(obj.x - obj.X);
K = obj.StiffnessMatrix;
% eval = eigs(K,reduced_dimension,'sm');
% for i = 1:reduced_dimension
%     if eval(i) < 0
%         disp(eval(i))
%         disp('eval of K negative')
%     end
% end



Mass = obj.M;
Mass = Mass(indLogical,indLogical);
K = K(indLogical,indLogical);

[V,D] = eigs(K,Mass,reduced_dimension,'sm');
[low_eig, permutation_indices] = sort(diag(D));
V = V(:,permutation_indices);

damping_eig = low_eig;
keep_ind = true(size(low_eig,1),1);
for i = 1:reduced_dimension
%     
%     current_def = obj.x - obj.X;
%     eig_mode = current_def;
%     eig_mode(indLogical) = current_def(indLogical) + V(:,i)*dt;
%     obj.SetCurrentState(eig_mode);
%     scal = 5;
%     offc = 0.05;
%     sky = surf(scal*(skyX-offc),scal*(skyY-offc),scal*(skyZ-offc),'LineStyle','none','FaceColor','flat','FaceLighting','gouraud'); %Make a sphere object.
%     
%     node = transpose(reshape(obj.x,3,[]));
%     tri1 = surftri(node,elem);
%     h = trimesh(tri1,obj.x(1:3:end),obj.x(2:3:end),obj.x(3:3:end));
%     %                 triplot(elems{i_filename},node(:,1),node(:,2),'Color',edge_colors{i_filename});
%     
%     axis equal
%     axis(axis_box)
%     % draw a horizontal reference line for droop test
%     %                 hl = refline(0,lowests{i_filename});
%     %                 hl.Color = edge_colors{i_filename};
%     
%     %                 ah = gca;
%     %                 phlist = findobj(ah,'Type','patch');
%     %                 ph = phlist(i_filename);
%     h.FaceLighting = 'flat';
%     h.FaceAlpha = 0.7;
%     h.FaceColor = 'interp';
%     %                 ph.EdgeLighting = 'gouraud';
%     %                 ph.BackFaceLighting = 'lit';
%     
%     %                 ph.EdgeColor = [0.9 0.9 1];
%     
%     %                 lh2 = camlight('headlight');
%     %                 lh.Position(2) = 3;
%     %                 lh.Position(3) = 3;
%     cmap = [0.8 0.9 1];
%     h.FaceVertexCData = repmat(cmap,N,1);
%     %                 colormap(autumn);
%     
%     shading faceted
%     %                 h.EdgeColor = edge_colors{i_filename};
%     %                 ph.EdgeLighting = 'none';
%     h.EdgeAlpha = 0.7;
%     %                 ph.LineWidth = 0.25;
%     %                 ph.AlignVertexCenters = 'on';
%     %                 whitebg(cmap)
%     
%     %                 origPos = [0,0,0]';
%     % %                 campos(origPos);
%     %                 origTarget = [0,0,-10000];
%     %                 camtarget(origTarget);
%     %                 camva(40)
%     ha.Visible = 'off';
%     ha.Clipping = 'off';
%     ha.Projection = 'perspective';
%     campos([ -2.1756    1.3321    2.1624]);
%     camtarget([0.0683    0.0419    0.0307])
%     %                 campos([0    2.0   2]);
%     %                 camtarget([0   -0.15    0.1])
%     camva(6.9295);
%     lighting gouraud
%     material metal
%     
%     colormap(ccmap)
%     hold on
%     %                 add_shadow()
%     
%     drawnow
%     
%     lh = camlight('headlight');
%     lh.Style = 'local';
%     lh.Position(3) = lh.Position(3) * 2;
%     ground = [];
%     color = [0.21 0.21 0.21];
%     nudge = 0.12;
%     background_color = [1 1 1];
%     fade = 'local';
%     T = [h]';
%     if isempty(ground)
%         minZ = inf;
%         for t = T'
%             V = t.Vertices;
%             minZ = min([V(:,3);minZ]);
%         end
%         ground = [0 0 -1 -nudge];
%     end
%     %             for i_filename = 1:length(trajectories)
%     %                 t = h{i_filename};
%     %                 t.EdgeColor = edge_colors{i_filename};
%     %                 t.FaceAlpha = 0.5;
%     %                 t.EdgeAlpha = 0.6;
%     %
%     %             end
%     h_shadow = [];
%     M = zeros([0,0,0]);
%     L = lh;
%     for t = T'
%         V = t.Vertices;
%         for l = L'
%             % plane equation
%             % 0 = ax + by + cz + d
%             light_pos = [l.Position strcmp(l.Style,'local')];
%             d = ground * light_pos';
%             shadow_mat = d*eye(4) - light_pos'*ground;
%             U = [V ones(size(V,1),1)]*shadow_mat';
%             U = bsxfun(@rdivide,U(:,1:3),U(:,4));
%             
%             hold on;
%             tsh = trisurf(t.Faces,U(:,1),U(:,2),U(:,3), ...
%                 'FaceColor',color, ...
%                 'DiffuseStrength',0,'SpecularStrength',0, ...
%                 'AmbientStrength',1, ...
%                 'EdgeColor','none');
%             %                     hold off;
%             switch fade
%                 case {'local','infinite'}
%                     D = matrixnormalize( ...
%                         sum(bsxfun(@times,tsh.Vertices(:,1:2),l.Position(1:2)),2));
%                     switch fade
%                         case 'infinite'
%                             D = 1.0-D;
%                     end
%                     tsh.FaceVertexCData = ...
%                         bsxfun(@plus,color,bsxfun(@times,D,background_color-color)/3-0.2);
%                     tsh.FaceColor = 'interp';
%             end
%             h_shadow = [h_shadow;tsh];
%             M(:,:,end+1) = shadow_mat;
%         end
%     end
%     print(fig,['Bar-oneside-fixed-coarse-step' num2str(dt) 'nonlinear-mode' num2str(i_mode)],'-dpng')
%     cla
%     hold on
%     obj.SetCurrentState(current_def);
    if low_eig(i) < 1e-3
        disp(low_eig(i))
%         low_eig(i) = 0;
        disp('eig below zero')
        keep_ind(i) = false;
    end
    
    if  ~isreal(low_eig(i))
        disp(low_eig(i))
%         low_eig(i) = 0;
        disp('eig imaginary')
        keep_ind(i) = false;
    end
    if damping_eig(i) < 0
        %         damping_eig(i) = -damping_eig(i);
        damping_eig(i) = 0;
    end
end
for i_eig = 1:length(obj.eig_ratios)
low_eig(i_eig) = low_eig(i_eig) * obj.eig_ratios(i_eig);
end
% keep_ind = diag(keep_ind);
V_pos = V(:,keep_ind);

reduced_dimension_pos = sum(keep_ind);

% for i = 1:reduced_dimension
% end

reduced_dq =  V_pos'*Mass*(obj.x(indLogical)-obj.X(indLogical));
reduced_v =  V_pos'*Mass*(v(indLogical));

B = -obj.a * Mass - obj.b * K;
reduced_B  = V_pos'*B*V_pos;

% for i = 1:reduced_dimension
%     if reduced_B(i,i) > 0
%         reduced_B(i,i) = 0;
%     end
% end


Eforce = obj.ElasticForce;
Eforce = Eforce(indLogical);

fExternal = Mass * obj.externalGravity(indLogical);


f = Eforce + fExternal + B*v(indLogical); % from column to row
% f = Eforce + fExternal;


rescaled_id = eye(obj.reduced_dimension);
rescaled_id(1:length(obj.eig_ratios),1:length(obj.eig_ratios)) = diag(obj.eig_ratios);

reduced_f =  rescaled_id *V_pos' * f;


% ERE

reduced_J = [sparse(reduced_dimension_pos,reduced_dimension_pos), speye(reduced_dimension_pos);...
    -diag(low_eig(keep_ind)), reduced_B];
% reduced_J = [sparse(reduced_dimension_pos,reduced_dimension_pos), speye(reduced_dimension_pos);...
%     -diag(low_eig(keep_ind)), -0*obj.b*diag(low_eig(keep_ind))];
reduced_du = [reduced_v; reduced_f];
reduced_g = reduced_du - reduced_J * [reduced_dq; reduced_v];
% eta = 2 ^ (-ceil(log2(norm(reduced_g,1))));
eta = 1;

reduced_J_tilde = sparse(size(reduced_J,1) + 1, size(reduced_J,1) + 1);
reduced_J_tilde(1:end-1,:) = [reduced_J, eta*reduced_g];
reduced_u_tilde = [[reduced_dq; reduced_v]; 1/eta];
X = expv(dt, reduced_J_tilde, reduced_u_tilde);

reduced_dq = X(1:(end-1)/2);
reduced_v = X((end-1)/2+1:end-1);

dq(indLogical) = V_pos * reduced_dq;
v(indLogical) = V_pos * reduced_v;
u_new = [dq; v];



end