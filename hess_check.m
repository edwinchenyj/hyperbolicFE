function hess
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


connectivity = [1:(n_total-1); 2:n_total]';
l = vecnorm(q(:,connectivity(:,1)) - q(:,connectivity(:,2)));

q = q(:);
v = zeros(size(q));

k = 10000000 * ones(size(l));


mass = 1;

springs = springs_particles(q, v, connectivity, k, 0, l,mass);

direction = rand(size(q));
direction = direction/vecnorm(direction);

epsilon = 1e-6;

springs_plus = springs_particles(q+direction*epsilon, v, connectivity, k, 0, l , mass);
springs_minus = springs_particles(q-direction*epsilon, v, connectivity, k, 0, l , mass);

force_plus = springs_plus.force;
force_minus = springs_minus.force;

grad_approx = (force_plus - force_minus)/2/epsilon;

grad_to_check = springs.force_derivative * direction;

display(max(abs(grad_approx - grad_to_check)))



function [x,y,z] = helix(n, r, h, d, p_per_turn)
theta = linspace(0,2*pi*n, n*p_per_turn);
x = r*cos(theta);
y = r*sin(theta);
z = linspace(0, h, n*p_per_turn);

end
end