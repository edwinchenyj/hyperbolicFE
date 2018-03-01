function [u_new, residual, step_length] = SemiImplicit_rescaled_fh( dt, u, obj, varargin)
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
obj.indLogical = indLogical;

Mass = obj.M;
Mass = Mass(indLogical,indLogical);

[Kh, Eforce] = obj.eigModification_nonlinear_subspace_fh;


fExternal = Mass * obj.externalGravity(indLogical);

% B = -obj.a * Mass - obj.b * K;


v = u(end/2 + 1:end);
f = Eforce -obj.a * Mass * v(indLogical) - obj.b * Kh(v(indLogical)) + fExternal;


%% version 1
% A = (Mass - dt * B + dt^2 * K);
% A = (Mass + dt^2 * K); % no rayleigh damping
rhs = dt * (f - dt * Kh(v(indLogical)));
dv_free = pcg(@A,rhs);


v_new_temp(indLogical) = v(indLogical) + dv_free;

residual = 1/2*(v_new_temp(indLogical) - v(indLogical) - dt * (obj.Minv*f))' * (v_new_temp(indLogical) - v(indLogical) - dt * (obj.Minv*f));
step_length = norm(dv_free);

v(indLogical) = v(indLogical) + dv_free;

% v(indLogical) = v(indLogical) + dv_free;
% max(v)

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

    function out = A(v)
        out = Mass*v - dt * (-obj.a * Mass*v - obj.b * Kh(v)) + (dt^2) * Kh(v);
    end
end