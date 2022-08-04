function der_sim


cam_pos = [5,5,1];        
    cam_el  = -.2;              
    cam_az  = -2;              
    
    cam_el_speed   = 0;           
    cam_az_speed   = 0;
    cam_fwd_speed  = 0;
    cam_side_speed = 0;
    
    unit_vec = @(el,az)[            
        cos(el)*cos(az), ...
        cos(el)*sin(az), ...
        sin(el)];
    
if(~ishandle(1))
fig = figure( ...
        'numbertitle','off', 'name','3d scene', ...
        'color',  [0.2    0.2    0.2], ...
        'closerequestfcn', @close_fcn);
else
    cla
    fig = gcf;
end

grid_dim = 50;
grid_spacing = 0.5;
grid_size = grid_dim/grid_spacing;
grid_tick = 0:grid_spacing:grid_dim;
    [ x_matrix, y_matrix ] = meshgrid(...
            grid_tick, ...
            grid_tick);
ax = axes;
ax.Projection = 'perspective';
ax.CameraViewAngle = 60;
ax.DataAspectRatio = [1,1,1];
ax.Visible = 'off';
plot_surf = surface(...
        'xdata', x_matrix, ...
        'ydata', y_matrix, ...
        'zdata', -1 * ones(size(x_matrix)), ...
        'cdata', 0.7*ones(grid_size,grid_size, 3), ...
        'facelighting', 'gouraud');
% set(plot_surf,'LineWidth',0.1);
set(plot_surf,'EdgeAlpha',0.3)

    % lights
%     light;
%     light('style', 'local');
lightangle(ax,-45,30)
% lightangle(ax,-30,30)
% lightangle(ax,-15,30)
% lighting gouraud
ax.CameraPosition = cam_pos;
        ax.CameraTarget = cam_pos + unit_vec(cam_el, cam_az);
n_loop = 1;
n_p_loop = 100;
n_total = n_loop * n_p_loop;
radius = 0.1;
height = 0.5;
[x,y,z] = helix(n_loop,radius,height,0.1, n_p_loop);
% x = linspace(0,1,n_total);
% y = zeros(size(x));
% z = -0.1*ones(size(x));

q = [x;y;-z];
translation = ax.CameraTarget' + [-0.5,-0.5,0.7]';
q = q + repmat(translation,[1,n_total]);


connectivity = [1:(n_total-1); 2:n_total]';
l = vecnorm(q(:,connectivity(:,1)) - q(:,connectivity(:,2)));
l = l;
q = q(:);
v = zeros(size(q));

k = 1e6 * ones(size(l));


mass = 1 * l(1);

springs = springs_particles(q, v, connectivity, k, 0, l,mass);
hold on

model = line(springs.q(1:3:end), springs.q(2:3:end), springs.q(3:3:end));
model.Color = [1 0.5 0];
model.LineWidth = 2;
model.Marker = 'none';
% model = plot3(1,1,1,'Color','red','LineWidth',2);
% 
%     model.XData = springs.q(1:3:end);
%     model.YData = springs.q(2:3:end);
%     model.ZData = springs.q(3:3:end);
dt = 0.001;
tsteps = 1000;

loop_on = true;
while loop_on
% for ti = 1:tsteps
    model.XData = springs.q(1:3:end);
    model.YData = springs.q(2:3:end);
    model.ZData = springs.q(3:3:end);
%     plot3(springs.q(1:3:end), springs.q(2:3:end), springs.q(3:3:end));
    drawnow
    indLogical = ones(size(springs.q));
    indLogical(1:3) = 0;
    indLogical = logical(indLogical);
    
    
    Mass = speye(nnz(indLogical));
    K = springs.stiffness_matrix;
    K = K(indLogical, indLogical);

%     B = -0.001 * Mass - 0.001 * K;

    gravity = zeros(size(indLogical));
    gravity(3:3:end) = -9.8;
    gravity = gravity(indLogical);

    force = springs.force;
    force = force(indLogical);
%% SI version 1
    A = (Mass + dt^2 * K);
    rhs = dt * ( force  + gravity - dt * K * springs.v(indLogical));
    dv = A\rhs;

    springs.v(indLogical) = springs.v(indLogical) + dv;

%%

    springs.q = springs.q + dt * springs.v;
    

end

function [x,y,z] = helix(n, r, h, d, p_per_turn)
theta = linspace(0,2*pi*n, n*p_per_turn);
x = r*cos(theta);
y = r*sin(theta);
z = linspace(0, h, n*p_per_turn);

end

    function close_fcn(~,~)
        loop_on = false;
        delete(fig);
    end

end