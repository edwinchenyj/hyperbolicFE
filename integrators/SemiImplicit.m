function [u_new, residual, step_length] = SemiImplicit( dt, u, obj, varargin)
% inputs:
%   dt: step size
%    u: current state
%  obj: the elastic obj
%  varargin: optional input for constraints
% notice sum(nnz(~constraint_indices)) = length(u)/2

if nargin == 3
    constraint_indices = false(size(u,1)/2,1);
elseif nargin == 4
    constraint_indices = varargin{1};
end


indLogical = ~constraint_indices;

K = obj.StiffnessMatrix;
Mass = obj.M;
Eforce = obj.ElasticForce;

Mass = Mass(indLogical,indLogical);
K = K(indLogical,indLogical);
B = -obj.a * Mass - obj.b * K;

Eforce = Eforce(indLogical);

fExternal = Mass * obj.externalGravity(indLogical);

vold = u(end/2 + 1:end);
v = vold;
f = Eforce + fExternal + B*v(indLogical);

%% version 1
A = (Mass - dt * B + dt^2 * K);
rhs = dt * (f - dt * K * v(indLogical));
dv_free = A\rhs;

% 
% v_new_temp(indLogical) = v(indLogical) + dv_free;
% 
% residual = 1/2*(v_new_temp(indLogical) - v(indLogical) - dt * (obj.Minv*f))' * (v_new_temp(indLogical) - v(indLogical) - dt * (obj.Minv*f));
% step_length = norm(dv_free);

v(indLogical) = v(indLogical) + dv_free;

%% version 2
% A = Mass\(Mass + dt^2 * K);
% rhs = v + dt * (Mass\f);
% v = A\rhs;

%% version 3
% A = Mass\(Mass + dt^2 * K);
% rhs = Mass\(Mass * v + dt * f);
% v = A\rhs;


u(end/2+1:end) = v;
u(1:end/2) = u(1:end/2) + dt * v;
obj.x = obj.X +  u(1:end/2);

obj.SetCurrentState(obj.x - obj.X);
u_new = u;


Eforce = obj.ElasticForce;
Eforce = Eforce(indLogical);
f = Eforce + fExternal + B*v(indLogical);
residual = norm((v(indLogical) - vold(indLogical) - dt * (Mass\f)));
step_length = norm(dv_free);

end