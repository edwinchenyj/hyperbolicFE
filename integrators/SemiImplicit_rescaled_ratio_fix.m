function u_new = SemiImplicit_rescaled_ratio_fix( dt, u, obj, varargin)
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

Mass = obj.M;
Mass = Mass(indLogical,indLogical);

[K, Eforce] = obj.eigModification_nonlinear_subspace_ratio_fix;

fExternal = Mass * obj.externalGravity(indLogical);



v = u(end/2 + 1:end);
f = Eforce + fExternal;


%% version 1
% A = (Mass - dt * B + dt^2 * K);
A = (Mass + dt^2 * K); % no rayleigh damping
rhs = dt * (f - dt * K * v(indLogical));
dv_free = A\rhs;

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
end