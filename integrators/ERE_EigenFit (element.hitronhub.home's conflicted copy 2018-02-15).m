function u_new = ERE_EigenFit( dt, u, obj, varargin)
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
Mass = obj.M;Mass = Mass(indLogical,indLogical);

K = K(obj.indLogical,obj.indLogical);
Mass = Mass(indLogical,indLogical);


[V_s, low_eig] = obj.EigenFit_subspace(K,Mass);

Eforce = obj.ElasticForce;
Eforce = Eforce(obj.indLogical);
Eforce =  Eforce + Mass * V_s * (diag(obj.eig_ratios) - eye(length(obj.eig_ratios))) * (V_s') * Eforce;
    ratio = spdiags(obj.eig_ratios,0,obj.eig_modes,obj.eig_modes);
    sI = speye(obj.eig_modes,obj.eig_modes);
    low_eig_diag = spdiags(low_eig,0,obj.eig_modes,obj.eig_modes);
rescaled_K = @(x) (K*x + (Mass * (V_s * (ratio - sI) * ((low_eig_diag) * (V_s') * Mass*x))));
    

fExternal = Mass * obj.externalGravity(indLogical);

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

    function out = J_tilde_h(in)
        len = length(in);
        x1 = in(1:(len-1)/2);
        x2 = in((len+1)/2:end-1);
        x_end = in(end);
    end
end