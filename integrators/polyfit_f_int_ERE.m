function u_new = polyfit_f_int_ERE( dt, u, obj, varargin)
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
K = obj.StiffnessMatrix;

Mass = obj.M;
Mass = Mass(indLogical,indLogical);
K = K(indLogical,indLogical);

Minv_K_new = polyvalm([obj.polyfit_p],Mass\K);

Eforce = obj.polyfit_force_approx;
% A = polyvalm([obj.polyfit_p(1:end-1)],Mass\K);
% Eforce = (A*Eforce(indLogical));
% Eforce = Mass*Minv_K_new*(K\Eforce(indLogical));
% [v_new, d_new] = eig(full(Minv_K_new),full(Mass));
K = Mass * Minv_K_new;


B = -obj.a * Mass - obj.b * K;



fExternal = Mass * obj.externalGravity(indLogical);


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

dq_old = u(1:end/2);
obj.polyfit_force_approx = obj.polyfit_force_approx - K * (dq(indLogical)-dq_old(indLogical));

end