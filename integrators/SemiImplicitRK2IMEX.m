function u_new = SemiImplicitRK2IMEX( dt, u, obj, varargin)
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


K = obj.DGElementStiffnessMatrix;
Mass = obj.M;
ElemEforce = obj.DGElementElasticForce;
Ki = obj.DGInterfaceStiffnessMatrix;
Mass = Mass(indLogical,indLogical);
K = K(indLogical,indLogical);
Ki = Ki(indLogical,indLogical);
B = -obj.a * Mass - obj.b * Ki;

ElemEforce = ElemEforce(indLogical);

InterfaceForce = obj.DGInterfaceElasticForce;
InterfaceForce = InterfaceForce(indLogical);

fExternal = Mass * obj.externalGravity(indLogical);

x = u(1:end/2);
v = u(end/2+1:end);
f = ElemEforce + InterfaceForce + fExternal + B*v(indLogical);

% rk2
v1 = v(indLogical) +  dt * (Mass\f);
x1 = x;
x1(indLogical) = x(indLogical) + dt * v1;
f1 = ElemEforce -Ki*x1 + fExternal + B*v(indLogical);



A = (Mass - dt * B + dt^2 * K);
f = 1/2 * (f + f1);
rhs = dt * (f - dt * K * v(indLogical));
dv_free = A\rhs;

v(indLogical) = v(indLogical) + dv_free;
u(end/2+1:end) = v;
u(1:end/2) = u(1:end/2) + dt * v;
obj.x = obj.X +  u(1:end/2);

obj.SetCurrentState(obj.x - obj.X);
u_new = u;
end