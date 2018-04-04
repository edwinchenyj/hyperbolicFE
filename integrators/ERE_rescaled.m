function [u_new, res, grad_size] = ERE_rescaled( dt, u, obj, varargin)
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

dq = u(1:end/2);
v = u(end/2+1:end);

N = obj.N;
nFixed = sum(constraint_indices)/Dim;

obj.x = obj.X + dq;

obj.SetCurrentState(obj.x - obj.X);

Mass = obj.M;
Mass = Mass(indLogical,indLogical);
[K, Eforce] = obj.eigModification_nonlinear_subspace;

fExternal = Mass * obj.externalGravity(indLogical);

B = -obj.a * Mass - obj.b * K;

f = Eforce + fExternal + B*v(indLogical); % from column to row

% ERE

J = [sparse(Dim*(N-nFixed),Dim*(N-nFixed)), speye(Dim*(N-nFixed)); -Mass\K, Mass\B];
du = [v(indLogical); Mass\f];
g = du - J * [dq(indLogical); v(indLogical)];
eta = 2 ^ (-ceil(log2(norm(g,1))));
% eta = 1;
J_tilde = sparse(size(J,1) + 1, size(J,1) + 1);
J_tilde(1:end-1,:) = [J, eta*g];
u_tilde = [[dq(indLogical); v(indLogical)]; 1/eta];
X = expv(dt, J_tilde, u_tilde);

dq(indLogical) = X(1:(end-1)/2);
v(indLogical) = X((end-1)/2+1:end-1);
u_new = [dq; v];

res =0;
grad_size = 0;

end