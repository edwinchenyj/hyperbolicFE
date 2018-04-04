function u_new = polyfit_BackwardEuler_modH_line_search( dt, u, obj, varargin)
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

it = 0;
dq = u(1:end/2);
v = u(end/2+1:end);
u_new = u;
v_new = u_new(end/2+1:end);

obj.x = obj.X + dq + dt * (v_new);

obj.SetCurrentState(obj.x - obj.X); % update deformation gradient
K_end = obj.StiffnessMatrix;
Mass = obj.M;
% Eforce_end = obj.ElasticForce;

Mass = Mass(indLogical,indLogical);
K_end = K_end(indLogical,indLogical);


Minv_K_new = polyvalm([obj.polyfit_p],Mass\K_end);
Eforce_end = obj.ElasticForce;
Eforce_end = Mass*Minv_K_new*(K_end\Eforce_end(indLogical));
K_end_modified = Mass * Minv_K_new;
% K_end_modified = obj.rest_stiffness;
% K_end_modified = K_end_modified(indLogical,indLogical);

B = -obj.a * Mass - obj.b * K_end;

Eforce_end = Eforce_end(indLogical);

fExternal = Mass * obj.externalGravity(indLogical);

f_end = Eforce_end + fExternal + B*v_new(indLogical);

N = obj.N;
nFixed = sum(constraint_indices)/obj.Dim;

% residual is the merit function in this case
residual_old = 1/2 *(dt * (Mass\f_end))' * (dt * (Mass\f_end));
Dv = -(speye(Dim*(N-nFixed)) + dt*dt*(Mass\K_end_modified) - dt*(Mass\B))\(v_new(indLogical) - v(indLogical) - dt * (Mass\f_end));


step_size = 1;
v_new_temp = v_new;
c1 = 1e-4;
c2 = 0.9;

v_new_temp(indLogical) = v_new(indLogical) + step_size*Dv;
dq_new_temp = dq + dt * (v_new_temp);
obj.SetCurrentState(dq_new_temp);
K_end = obj.StiffnessMatrix;
Mass = obj.M;

Mass = Mass(indLogical,indLogical);
K_end = K_end(indLogical,indLogical);
Minv_K_new = polyvalm([obj.polyfit_p],Mass\K_end);
Eforce_end = obj.ElasticForce;
Eforce_end = Mass*Minv_K_new*(K_end\Eforce_end(indLogical));
Eforce_end = Eforce_end(indLogical);
f_end = Eforce_end + fExternal + B*v_new_temp(indLogical);


% merit function
residual = 1/2*(v_new_temp(indLogical) - v(indLogical) - dt * (Mass\f_end))' * (v_new_temp(indLogical) - v(indLogical) - dt * (Mass\f_end));

% wolfe conditions
% while (residual > residual_old + c1 * step_size * residual * 2) || (residual < c2 * residual_old)
while (residual > residual_old + c1 * step_size * residual * 2)
    residual_old = residual;
    it = it + 1;
    if it == 40
        error('line search more than 40 iterations.');
    end
    step_size = 0.9 * step_size;
    v_new_temp(indLogical) = v_new(indLogical) + step_size*Dv;
dq_new_temp = dq + dt * (v_new_temp);
obj.SetCurrentState(dq_new_temp);
K_end = obj.StiffnessMatrix;
Mass = obj.M;

Mass = Mass(indLogical,indLogical);
K_end = K_end(indLogical,indLogical);
Minv_K_new = polyvalm([obj.polyfit_p],Mass\K_end);
Eforce_end = obj.ElasticForce;
Eforce_end = Mass*Minv_K_new*(K_end\Eforce_end(indLogical));
Eforce_end = Eforce_end(indLogical);
f_end = Eforce_end + fExternal + B*v_new_temp(indLogical);

    % merit function
    residual = 1/2*(v_new_temp(indLogical) - v(indLogical) - dt * (Mass\f_end))' * (v_new_temp(indLogical) - v(indLogical) - dt * (Mass\f_end));
    
end

v_new = v_new_temp;

u_new(1:end/2) = dq + dt * (v_new);
u_new(end/2+1:end) = v_new;


end