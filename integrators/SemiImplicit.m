function u_new = SemiImplicit( dt, u, obj, varargin)
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

fExternal = Mass * obj.externalGravity;


f = Eforce + fExternal + B*u(end/2 + 1:end);

A = (Mass - dt * B + dt^2 * K);
v_free = u(end/2 + 1:end);
rhs = dt * (f - dt * K * v_free);
dv_free = A\rhs;

u(end/2 +1 :end) = u(end/2 +1 :end) + dv_free;
u(1:end/2) = u(1:end/2) + dt * u(end/2+1:end);
obj.x(indLogical) = obj.X(indLogical) +  u(1:end/2);

obj.SetCurrentState(obj.x - obj.X);

u_new = u;
end