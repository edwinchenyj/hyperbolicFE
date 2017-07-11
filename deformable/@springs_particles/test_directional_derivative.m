function test_directional_derivative
    clear all
close all
clc
% initialize n particles on a vertial chain with seperation l(i)
n = 10;
q = zeros( 3 * n, 1 );


% fixing a random seed
rng(1)
% initialize m springs
m = n - 1;
% stiffness
k = 1*ones(m,1);
% linear damping constants
k_damping = 0 * rand(m,1);
springs = zeros(m,2);
% each row contains the indices of the particles connected by that spring
springs(:,1) = 1:(n-1); springs(:,2) = 2:n; 
% original length of each spring
l = rand(n-1,1);

% mass of each particle
mass = ones(3*n,1);
% mass = rand(n,1);
% mass = repmat(mass,1,3);
% mass = mass(:);
% initial velocity of each particle
v = rand(3*n,1) * .1;
v(1:3) = zeros(3,1);
% initial positions
for i = 2:n
    q( 3 * i - 1) = q( 3 * (i-1) - 1 ) - l(i-1);
end

init_state = [q;v];
% constructor 
% obj = springs_particles(q,v,springs,k,k_damping,l,mass)
SPobject = springs_particles(init_state(1:3*n),init_state(3*n+1:end),springs,k,k_damping,l,mass);

for i = 1 : 8
dx = normc(rand(3*n,1));
f0 = SPobject.force();
K = SPobject.stiffness_matrix();
h = 10^(-i);
SPobject.update(q + h * dx, v);
f1 = SPobject.force();

(f1 - f0)/h
K * dx
norm((f1 - f0)/h)
norm(K * dx)
end
end
