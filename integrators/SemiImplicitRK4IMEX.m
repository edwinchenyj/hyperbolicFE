function u_new = SemiImplicitRK4IMEX( dt, u, obj, varargin)
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

% rk4
v1 = v(indLogical) + 1/2 * dt * (Mass\f);
x1 = x;
x1(indLogical) = x(indLogical) + 1/2 * dt * v;
f1 = ElemEforce -Ki*x1 + fExternal + B*v(indLogical);
v2 = v(indLogical) + 1/2 * dt * (Mass\f1);
x2 = x;
x2(indLogical) = x(indLogical) + 1/2 * dt * v2;
f2 = ElemEforce -Ki*x2 + fExternal + B*v2(indLogical);
v3 = v(indLogical) + dt * (Mass\f2);
x3 = x;
x3(indLogical) = x(indLogical) + dt * v3;
f3 = ElemEforce -Ki*x3 + fExternal + B*v3(indLogical);

A = (Mass - dt * B + dt^2 * K);
f = ElemEforce - Ki *( 1/6 * x(indLogical) + 1/3 * x1(indLogical) + 1/3 * x2(indLogical) + 1/6 * x3(indLogical)) + fExternal...
    + B*( 1/6 * v(indLogical) + 1/3 * v1(indLogical) + 1/3 * v2(indLogical) + 1/6 * v3(indLogical));
rhs = dt * (f - dt * K * v(indLogical));
dv_free = A\rhs;

v(indLogical) = v(indLogical) + dv_free;
u(end/2+1:end) = v;
u(1:end/2) = u(1:end/2) + dt * v;
obj.x = obj.X +  u(1:end/2);

obj.SetCurrentState(obj.x - obj.X);
u_new = u;
end