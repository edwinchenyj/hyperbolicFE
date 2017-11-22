function u_new = SemiImplicitModified( dt, u, obj, varargin)
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

v = u(end/2 + 1:end);
f = Eforce + fExternal + B*v(indLogical);

A = (Mass - dt * B + 0.5*dt^2 * K); % todo: update B
rhs = dt * (f - 0.5*dt * K * v(indLogical));
dv_free = A\rhs;

u(1:end/2) = u(1:end/2) + 0.5*dt * v;
v(indLogical) = v(indLogical) + dv_free;
u(end/2+1:end) = v;
u(1:end/2) = u(1:end/2) + 0.5*dt * v;
obj.x = obj.X +  u(1:end/2);

obj.SetCurrentState(obj.x - obj.X);
u_new = u;
end