function u_new = reduced_extended_ERE( dt, u, obj, varargin)
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
reduced_dimension = obj.reduced_dimension;


dq = u(1:end/2);
v = u(end/2+1:end);

N = obj.N;
nFixed = sum(constraint_indices)/Dim;

obj.x = obj.X + dq;

obj.SetCurrentState(obj.x - obj.X);
K = obj.StiffnessMatrix;
% eval = eigs(K,reduced_dimension,'sm');
% for i = 1:reduced_dimension
%     if eval(i) < 0
%         disp(eval(i))
%         disp('eval of K negative')
%     end
% end



Mass = obj.M;
Mass = Mass(indLogical,indLogical);
K = K(indLogical,indLogical);

[V,D] = eigs(K,Mass,reduced_dimension,'sm');
[low_eig, permutation_indices] = sort(diag(D));
V = V(:,permutation_indices);

damping_eig = low_eig;
keep_ind = true(size(low_eig,1),1);
for i = 1:reduced_dimension

    if low_eig(i) < 1e-3
        disp(low_eig(i))
%         low_eig(i) = 0;
        disp('eig below zero')
        keep_ind(i) = false;
    end
    
    if  ~isreal(low_eig(i))
        disp(low_eig(i))
%         low_eig(i) = 0;
        disp('eig imaginary')
        keep_ind(i) = false;
    end
    if damping_eig(i) < 0
        %         damping_eig(i) = -damping_eig(i);
        damping_eig(i) = 0;
    end
end
% keep_ind = diag(keep_ind);
V_pos = V(:,keep_ind);

reduced_dimension_pos = sum(keep_ind);

compressed_dq =  V_pos * (V_pos')*Mass*(obj.x(indLogical)-obj.X(indLogical));
compressed_v =  V_pos *(V_pos')*Mass*(v(indLogical));

B = -obj.a * Mass - obj.b * K;
compressed_Minv_B  = V_pos *(V_pos')*B*V_pos * (V_pos')* Mass;
compressed_Minv_K = V_pos *(V_pos')*K*V_pos * (V_pos') * Mass;


% for i = 1:reduced_dimension
%     if reduced_B(i,i) > 0
%         reduced_B(i,i) = 0;
%     end
% end

% % dq(indLogical) = compressed_dq;
% obj.SetCurrentState(dq);
Eforce = obj.ElasticForce;
Eforce = Eforce(indLogical);

fExternal = Mass * obj.externalGravity(indLogical);


f = Eforce + fExternal + B*v(indLogical); 
% f = Eforce + fExternal;

compressed_Minv_f =  V_pos * (V_pos') * f;

% ERE

J = [sparse(Dim*(N-nFixed),Dim*(N-nFixed)), speye(Dim*(N-nFixed)); -compressed_Minv_K, compressed_Minv_B];
% J = [sparse(Dim*(N-nFixed),Dim*(N-nFixed)), speye(Dim*(N-nFixed)); -Mass\K, Mass\B];
du = [compressed_v; compressed_Minv_f];
% du = [compressed_v; Mass\f];
% du = [v(indLogical); compressed_Minv_f];
g = du - J * [dq(indLogical); v(indLogical)];
g = du - J * [compressed_dq; compressed_v];
eta = 2 ^ (-ceil(log2(norm(g,1))));
% eta = 1;
J_tilde = sparse(size(J,1) + 1, size(J,1) + 1);
J_tilde(1:end-1,:) = [J, eta*g];
% u_tilde = [[dq(indLogical); v(indLogical)]; 1/eta];
u_tilde = [[compressed_dq; compressed_v]; 1/eta];
X = expv(dt, J_tilde, u_tilde);

dq(indLogical) = X(1:(end-1)/2);
v(indLogical) = X((end-1)/2+1:end-1);
u_new = [dq; v];


end