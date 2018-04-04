function u_new = ERE_rescaled_ratio_fixed( dt, u, obj, varargin)
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
[K, Eforce] = obj.eigModification_nonlinear_subspace_ratio_fix;

fExternal = Mass * obj.externalGravity(indLogical);


% 
% [V,D] = eigs(K,Mass,reduced_dimension,'sm');
% [low_eig, permutation_indices] = sort(diag(D));
% V = V(:,permutation_indices);
% 
% damping_eig = low_eig;
% keep_ind = true(size(low_eig,1),1);
% for i = 1:reduced_dimension
% 
%     if low_eig(i) < 1e-3
%         disp(low_eig(i))
% %         low_eig(i) = 0;
%         disp('eig below zero')
%         keep_ind(i) = false;
%     end
%     
%     if  ~isreal(low_eig(i))
%         disp(low_eig(i))
% %         low_eig(i) = 0;
%         disp('eig imaginary')
%         keep_ind(i) = false;
%     end
%     if damping_eig(i) < 0
%         %         damping_eig(i) = -damping_eig(i);
%         damping_eig(i) = 0;
%     end
% end
% % keep_ind = diag(keep_ind);
% V_pos = V(:,keep_ind);
% 
% reduced_dimension_pos = sum(keep_ind);
% 
% compressed_dq =  V_pos * (V_pos')*Mass*(obj.x(indLogical)-obj.X(indLogical));
% compressed_v =  V_pos *(V_pos')*Mass*(v(indLogical));
% 
% B = -obj.a * Mass - obj.b * K;
% compressed_Minv_B  = V_pos *(V_pos')*B*V_pos * (V_pos')* Mass;
% compressed_Minv_K = V_pos *(V_pos')*K*V_pos * (V_pos') * Mass;
% 
% 
% % for i = 1:reduced_dimension
% %     if reduced_B(i,i) > 0
% %         reduced_B(i,i) = 0;
% %     end
% % end
% 
% % % dq(indLogical) = compressed_dq;
% % obj.SetCurrentState(dq);
% Eforce = obj.ElasticForce;
% Eforce = Eforce(indLogical);
% 
% fExternal = Mass * obj.externalGravity(indLogical);
% 

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

end