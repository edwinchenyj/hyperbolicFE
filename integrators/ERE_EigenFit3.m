function [u_new, res, gradient_size] = ERE_EigenFit2( dt, u, obj, varargin)
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
len = sum(indLogical);
Dim = obj.Dim;

dq = u(1:end/2);
v = u(end/2+1:end);

N = obj.N;
nFixed = sum(constraint_indices)/Dim;

obj.x = obj.X + dq;

obj.SetCurrentState(obj.x - obj.X);

K = obj.StiffnessMatrix;
Mass = obj.M;

K = K(obj.indLogical,obj.indLogical);
Mass = Mass(indLogical,indLogical);


[V_s, low_eig] = obj.EigenFit_subspace(K,Mass);

Eforce = obj.ElasticForce;
Eforce = Eforce(obj.indLogical);
Eforce =  Eforce + Mass * V_s * (diag(obj.eig_ratios) - eye(length(obj.eig_ratios))) * (V_s') * Eforce;
    ratio = spdiags(obj.eig_ratios,0,obj.eig_modes,obj.eig_modes);
    sI = speye(obj.eig_modes,obj.eig_modes);
    low_eig_diag = spdiags(low_eig,0,obj.eig_modes,obj.eig_modes);
    rescaled_const = (ratio - sI) * low_eig_diag;
    rescaled_V_s = rescaled_const * (V_s');
% rescaled_K = @(x) K*x + Mass * (V_s * (rescaled_V_s * (Mass*x)));

fExternal = Mass * obj.externalGravity(indLogical);


% [K, Eforce] = obj.eigModification_nonlinear_subspace;
% 
% fExternal = Mass * obj.externalGravity(indLogical);
% 
% B = -obj.a * Mass - obj.b * K;
% 
f = Eforce + fExternal + B(v(indLogical));

% ERE

% J = [sparse(Dim*(N-nFixed),Dim*(N-nFixed)), speye(Dim*(N-nFixed)); -Mass\K, Mass\B];
du = [v(indLogical); Mass\f];
% g = du - J * [dq(indLogical); v(indLogical)];
% g = du - [v(indLogical);-Mass\rescaled_K(dq(indLogical))+Mass\B(v(indLogical))];
g = du - [v(indLogical);-Mass\rescaled_K(dq(indLogical))];

eta = 2 ^ (-ceil(log2(norm(g,1))));
% eta = 1;
% J_tilde = sparse(size(J,1) + 1, size(J,1) + 1);
% J_tilde(1:end-1,:) = [J, eta*g];
u_tilde = [[dq(indLogical); v(indLogical)]; 1/eta];
% norm_approx = norm(K,'inf');
norm_approx = max(J_tilde(ones(length(g)+1,1)));
X = expv_fh(dt, @J_tilde, u_tilde,norm_approx);

dq(indLogical) = X(1:(end-1)/2);
v(indLogical) = X((end-1)/2+1:end-1);
u_new = [dq; v];

res = 0;
gradient_size = 0;

    function output = J(input)
        x1 = input(1:end/2);
        x2 = input(end/2+1:end);
        output = zeros(size(input));
        output(1:end/2) = x2;
        output(end/2+1:end) = -Mass\rescaled_K(x1) + Mass\B(x2);
    end

    function output = J_tilde(input)
%         x1 = input(1:(end-1)/2);
%         x2 = input((end+1)/2:end-1);
        xend = input(end);
        output = zeros(size(input));
        output(1:(end-1)) = J(input(1:end-1)) + eta*g*xend;
%         output((end+1)/2:end-1) = -Mass\rescaled_K(x1)+eta*g(end/2+1:end)*xend;
        
    end

    function output = B(input)
        output = -obj.a * Mass * input - obj.b * rescaled_K(input);
    end

    function output = rescaled_K(input)
        term1 = K*input;
        term2 = Mass*input;
        term2 = rescaled_V_s* term2;
        term2 = V_s * term2;
        term2 = Mass * term2;
        
        output = term1 + term2;
%     K*x + Mass * (V_s * (rescaled_V_s * (Mass*x)))
    end
end