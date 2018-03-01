function name = simulation3D_motions_bar_mode_visualization(varargin)
close all
clf

draw = true;
rerun_flag = true;
save_state = true;

dt = 1/100;
T = 4;
tsteps = 100*T;

fs = filesep;

mesh_shape = 'small_bar';
h = 0.2;
simulation_type = 'CG';



solver = 'SI';

Y = 6000; % Young's modululs
P = 0.45; % Poisson ratio
rho = 1000; % density
a = 0.0; % rayleigh damping
b = 0.00;
material_type = 'neo-hookean'; % choice: 'linear', 'neo-hookean'

axis_box = [-1 1 -1 1 -1 1]/4;


gravity = 'off';
mode = 'n';
constraint = 'n';
deformation_scale = 'n';
DGeta = 'n';

% parse input
i_arg = 1;
while (i_arg <= nargin)
    switch varargin{i_arg}
        case 'draw'
            i_arg = i_arg + 1;
            draw = varargin{i_arg};
        case 'rerun_flag'
            i_arg = i_arg + 1;
            rerun_flag = varargin{i_arg};
        case 'save_state'
            i_arg = i_arg + 1;
            save_state = varargin{i_arg};
        case 'dt'
            i_arg = i_arg + 1;
            dt = varargin{i_arg};
        case 'T'
            i_arg = i_arg + 1;
            T = varargin{i_arg};
        case 'mesh_shape'
            i_arg = i_arg + 1;
            mesh_shape = varargin{i_arg};
        case 'h'
            i_arg = i_arg + 1;
            h = varargin{i_arg};
        case 'simulation_type'
            i_arg = i_arg + 1;
            simulation_type = varargin{i_arg};
        case 'DGeta'
            i_arg = i_arg + 1;
            DGeta = varargin{i_arg};
        case 'solver'
            i_arg = i_arg + 1;
            solver = varargin{i_arg};
        case 'Y'
            i_arg = i_arg + 1;
            Y = varargin{i_arg};
        case 'P'
            i_arg = i_arg + 1;
            P = varargin{i_arg};
        case 'rho'
            i_arg = i_arg + 1;
            rho = varargin{i_arg};
        case 'a'
            i_arg = i_arg + 1;
            a = varargin{i_arg};
        case 'b'
            i_arg = i_arg + 1;
            b = varargin{i_arg};
        case 'material'
            i_arg = i_arg + 1;
            material_type = varargin{i_arg};
        case 'axis_box'
            i_arg = i_arg + 1;
            axis_box = varargin{i_arg};
        case 'constraint'
            i_arg = i_arg + 1;
            constraint = varargin{i_arg};
        case 'mode'
            i_arg = i_arg + 1;
            mode = varargin{i_arg};
        case 'gravity'
            i_arg = i_arg + 1;
            gravity = varargin{i_arg};
        case 'deformation_scale'
            i_arg = i_arg + 1;
            deformation_scale = varargin{i_arg};
            
    end
    i_arg = i_arg + 1;
end

tsteps = T/dt;

%% process the high res mesh first
h_original = h;
h = 1/10;

meshname = sprintf('mesh_data%c%s_h_%.2d',fs,mesh_shape, h);

if exist([meshname '.mat'], 'file') ~= 2
    error('mesh does not exist')
    
else
    load([meshname '.mat'], 'nodeM', 'elem');
    
end
% nodeM = nodeM(:,[2 1]);
% elem = elem(:,[1 3 2]);
N = size(nodeM,1);
% 
% loaddir = sprintf('sim_data%csimulation3D_rotate_bar_%s_%s', fs,material_type, mesh_shape);
% simdir = strcat(loaddir,fs,'SI','_','CG',...
%     '_constraint_',constraint,...
%     '_h',num2str(h),...
%     '_Y',num2str(1e5),...
%     '_P',num2str(0.45),...
%     '_rho',num2str(1000),...
%     '_a',num2str(0.01),...
%     '_b',num2str(0.00),...
%     '_dt',num2str(dt),...
%     '_def-scl',num2str(deformation_scale));
% loaded =load([simdir fs 'trajectory.mat'],'obj');
% 
% dirname = sprintf('sim_data%c%s_%s_%s', fs, mfilename,material_type, mesh_shape);
% 
% % fix nFixed number of highest point

[~, ind] = sortrows(nodeM,1);

% put the nodes in order (in z)
nodeM = nodeM(ind,:);

% create the rank for the original z
[~,R] = sort(ind);
% the new element list
elem = R(elem(:,:));

% nodeM = loaded.obj.nodeM;
% elem = loaded.obj.elem;

left_points = find(abs(nodeM(:,1)-0) < 0.01);
right_points = find(abs(nodeM(:,1)-1/5) < 0.01/5);

% get logical indices for the left and right nodes

indLeft = left_points;
indRight = right_points;

nFixed = length(indRight);

N = size(nodeM,1); % number of nodes
indAll = 1:N;
indRemove = indAll([indRight(:)]);
indLogical = logical(ones(3,N));
indLogical(:,indRemove) = logical(0);
indLogical = indLogical(:);

indLeftLogical = logical(zeros(3,N));
indLeftLogical(:,indLeft) = logical(1);
indLeftLogical = indLeftLogical(:);

indRightLogical = logical(zeros(3,N));
indRightLogical(:,indRight) = logical(1);
indRightLogical(:);


obj = elasticTetObj(nodeM, elem);
switch material_type
    case 'linear'
        obj.SetMaterial( Y, P, rho, 2, a, b); % set the tri to linear
        %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
    case 'neo-hookean'
        obj.SetMaterial( Y, P, rho, 1, a, b); % set the tri to neo-hookean
end


Dx = 0*rand(obj.Dim*N,1); % displacement field. set zero for rest state
obj.SetCurrentState(Dx);
%
M = obj.M;
K = obj.StiffnessMatrix;
M = 1/2 * (M + M');
K = 1/2 * (K + K');
M = M(indLogical,indLogical);
K = K(indLogical,indLogical);
%         K = K(~ind_fix,~ind_fix); % extract the non-fixed part

[V,D] = eigs(K,M,6,'sm');
[low_eig, permutation_indices] = sort(diag(D));
eigv = V(:,permutation_indices);
high_res_low_eig = low_eig;


v = zeros(size(Dx));
ha = obj.init_vis;
obj.simple_vis(obj.vis_handle);

xlim = axis_box(1:2);
ylim = axis_box(3:4);
zlim = axis_box(5:6);
u = [Dx; v];
obj.vis_handle.XLim = xlim;
obj.vis_handle.YLim = ylim;
obj.vis_handle.ZLim = zlim;
ha = obj.vis_handle;
camP = campos(gca);
camV = camva(gca);

fig = gcf;

fig.Children.Visible = 'off';
fig.Children.Clipping = 'off';
fig.Children.Projection = 'perspective';
fig.Position(3:4) = fig.Position(3:4) *2;
fig.Position(1:2) = [1300,50];
                ccmap = gray - 0.8;
                ind = ccmap < 0;
                ccmap(ind) = 0;
skymap = colormap(ccmap); %We'll use the colormap "winter" which is defined by matlab.

% Here's a gradient I made up that's a bit more horizon-like. You can
% experiment.
% flipud([logspace(-0.5,0,1000)',logspace(-0.5,0,1000)',logspace(-0.001,0,1000)']);

[skyX,skyY,skyZ] = sphere(200); %create the surface data for a sphere.
cla
hold on
% rate to draw the scene
sim_rate = round(1/dt);
patch_color = [0.8,0.9,0.5];

for i_mode = 1:6
    Dx(indLogical) = eigv(:,i_mode)/20;
    obj.SetCurrentState(Dx);
                scal = 5;
            offc = 0.05;
            sky = surf(scal*(skyX-offc),scal*(skyY-offc),scal*(skyZ-offc),'LineStyle','none','FaceColor','flat','FaceLighting','gouraud'); %Make a sphere object.
            
 node = transpose(reshape(obj.x,3,[]));
                tri1 = surftri(node,elem);
                h = trimesh(tri1,obj.x(1:3:end),obj.x(2:3:end),obj.x(3:3:end));
%                 triplot(elems{i_filename},node(:,1),node(:,2),'Color',edge_colors{i_filename});
                
                axis equal
                axis(axis_box)
                % draw a horizontal reference line for droop test
                %                 hl = refline(0,lowests{i_filename});
                %                 hl.Color = edge_colors{i_filename};
                
%                 ah = gca;
%                 phlist = findobj(ah,'Type','patch');
%                 ph = phlist(i_filename);
                h.FaceLighting = 'flat';
                h.FaceAlpha = 0.7;
                h.FaceColor = 'interp';
                %                 ph.EdgeLighting = 'gouraud';
                %                 ph.BackFaceLighting = 'lit';

                %                 ph.EdgeColor = [0.9 0.9 1];

                %                 lh2 = camlight('headlight');
                %                 lh.Position(2) = 3;
                %                 lh.Position(3) = 3;
                cmap = [0.8 0.9 1];
                h.FaceVertexCData = repmat(cmap,N,1);
                %                 colormap(autumn);
                
                shading faceted
%                 h.EdgeColor = edge_colors{i_filename};
%                 ph.EdgeLighting = 'none';
                h.EdgeAlpha = 0.7;
                %                 ph.LineWidth = 0.25;
                %                 ph.AlignVertexCenters = 'on';
%                 whitebg(cmap)
                
                %                 origPos = [0,0,0]';
                % %                 campos(origPos);
                %                 origTarget = [0,0,-10000];
                %                 camtarget(origTarget);
                %                 camva(40)
                ha.Visible = 'off';
                ha.Clipping = 'off';
                ha.Projection = 'perspective';
                campos([   -0.8334    1.5512    0.6300]);
                camtarget([0.1045    0.0087   -0.0369])
%                 campos([0    2.0   2]);
%                 camtarget([0   -0.15    0.1])
                camva(6.9295);
                lighting gouraud
                material metal

                colormap(ccmap)
                hold on
%                 add_shadow()

                drawnow
    
            lh = camlight('headlight');
            lh.Style = 'local';
            lh.Position(3) = lh.Position(3) * 2;
            ground = [];
            color = [0.21 0.21 0.21];
            nudge = 0.12;
            background_color = [1 1 1];
            fade = 'local';
            T = [h]';
            if isempty(ground)
                minZ = inf;
                for t = T'
                    V = t.Vertices;
                    minZ = min([V(:,3);minZ]);
                end
                ground = [0 0 -1 -nudge];
            end
%             for i_filename = 1:length(trajectories)
%                 t = h{i_filename};
%                 t.EdgeColor = edge_colors{i_filename};
%                 t.FaceAlpha = 0.5;
%                 t.EdgeAlpha = 0.6;
%             
%             end
            h_shadow = [];
            M = zeros([0,0,0]);
            L = lh;
            for t = T'
                V = t.Vertices;
                for l = L'
                    % plane equation
                    % 0 = ax + by + cz + d
                    light_pos = [l.Position strcmp(l.Style,'local')];
                    d = ground * light_pos';
                    shadow_mat = d*eye(4) - light_pos'*ground;
                    U = [V ones(size(V,1),1)]*shadow_mat';
                    U = bsxfun(@rdivide,U(:,1:3),U(:,4));
                    
                    hold on;
                    tsh = trisurf(t.Faces,U(:,1),U(:,2),U(:,3), ...
                        'FaceColor',color, ...
                        'DiffuseStrength',0,'SpecularStrength',0, ...
                        'AmbientStrength',1, ...
                        'EdgeColor','none');
%                     hold off;
                    switch fade
                        case {'local','infinite'}
                            D = matrixnormalize( ...
                                sum(bsxfun(@times,tsh.Vertices(:,1:2),l.Position(1:2)),2));
                            switch fade
                                case 'infinite'
                                    D = 1.0-D;
                            end
                            tsh.FaceVertexCData = ...
                                bsxfun(@plus,color,bsxfun(@times,D,background_color-color)/3-0.2);
                            tsh.FaceColor = 'interp';
                    end
                    h_shadow = [h_shadow;tsh];
                    M(:,:,end+1) = shadow_mat;
                end
            end
            print(fig,['Bar-oneside-fixed-coarse-Y' num2str(Y) '-mode' num2str(i_mode)],'-dpng')
            cla
            hold on
end

end