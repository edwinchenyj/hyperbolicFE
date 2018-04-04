function [u_new, residual, step_length] = SemiImplicit_EigenFit( dt, u, obj, varargin)
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

K = obj.StiffnessMatrix;
Mass = obj.M;

K = K(obj.indLogical,obj.indLogical);
Mass = Mass(indLogical,indLogical);

[V_s, low_eig] = obj.EigenFit_subspace(K,Mass);

Eforce = obj.ElasticForce;
Eforce = Eforce(obj.indLogical);
% Eforce =  Eforce + Mass * V_s * (diag(obj.eig_ratios) - eye(length(obj.eig_ratios))) * (V_s') * Eforce;
    ratio = spdiags(obj.eig_ratios,0,obj.eig_modes,obj.eig_modes);
    sI = speye(obj.eig_modes,obj.eig_modes);
    
    low_eig_diag = spdiags(low_eig,0,obj.eig_modes,obj.eig_modes);
    rescaled_const = (ratio - sI) * low_eig_diag;
    rescaled_V_s = rescaled_const * (V_s');
rescaled_K = @(x) K*x + Mass * (V_s * (rescaled_V_s * (Mass*x)));
    

fExternal = Mass * obj.externalGravity(indLogical);

% B = -obj.a * Mass - obj.b * K;


vold = u(end/2 + 1:end);
v = vold;
f = Eforce -obj.a * Mass * v(indLogical) - obj.b * rescaled_K(v(indLogical)) + fExternal;


%% version 1
% A = (Mass - dt * B + dt^2 * K);
% A = (Mass + dt^2 * K); % no rayleigh damping
rhs = dt * (f - dt * rescaled_K(v(indLogical)));
% dv_free = A\rhs;
tol = 1e-4;  
maxit = 10000;  
% M1 = spdiags(ones(sum(indLogical),1),0,sum(indLogical),sum(indLogical));
% tic
[dv_free, flag] = pcg(@A,rhs,tol,maxit);
% toc
% v_new_temp(indLogical) = v(indLogical) + dv_free;
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
% 
% Eforce = obj.eigModification_nonlinear_subspace_f_only;
% 
% f = Eforce -obj.a * Mass * v(indLogical) - obj.b * rescaled_K(v(indLogical)) + fExternal;

residual = norm((v(indLogical) - vold(indLogical) - dt * (Mass\f)));
step_length = norm(dv_free);



    function out = A(x)
        out = Mass*x - dt * (-obj.a * Mass*x) + ((dt^2) + dt * obj.b) * rescaled_K(x); 
    end
end