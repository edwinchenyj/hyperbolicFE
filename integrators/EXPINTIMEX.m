function u_new = EXPINTIMEX( dt, u, obj, varargin)
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

x = u(1:end/2);
v = u(end/2+1:end);
x_free = x(indLogical);
v_free = v(indLogical);

Ke = obj.DGElementStiffnessMatrix;
Mass = obj.M;
% Eforce = obj.ElasticForce;

Mass = Mass(indLogical,indLogical);
Ke = Ke(indLogical,indLogical);

Ki = obj.DGInterfaceStiffnessMatrix;
Ki = Ki(indLogical,indLogical);
B = -obj.a * Mass - obj.b * Ki;

A = sparse(size(Ki,1)*2,size(Ki,1)*2);
A(end/2+1:end,1:end/2) = -Mass\Ki;
A(end/2+1:end,end/2+1:end) = Mass\B;

if isempty(obj.InterfaceExp)
    obj.InterfaceExp = expm(dt*A);
end
X = obj.InterfaceExp * [x_free;v_free];

fe = obj.DGElementElasticForce;
fe = fe(indLogical);

fExternal = Mass * obj.externalGravity(indLogical);

rhs = X(end/2+1:end) + dt *(Mass\(fe+fExternal));
LHS = speye(size(x_free,1));
% LHS(1:end/2,end/2+1:end) = -dt * speye(size(x_free,1),size(x_free,1));
LHS = LHS  + dt*dt * (Mass\Ke);

v_free_new = LHS\rhs;

% x(indLogical) = u_free_new(1:end/2);
v(indLogical) = v_free_new;


u(end/2+1:end) = v;
u(1:end/2) = u(1:end/2) + dt * v;
obj.x = obj.X +  u(1:end/2);

% u = [x;v];

% obj.x = obj.X +  u(1:end/2);

obj.SetCurrentState(obj.x - obj.X);
u_new = u;
end